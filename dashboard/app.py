"""
app.py — FastAPI backend for the reaction-route dashboard, VIEW (read-only) build.

This is the public/viewer counterpart of scripts/dashboard/to_edit/app.py. It is
deliberately self-contained (its own copy of chem.py) so the whole to_view/ folder
can be lifted out and deployed as a standalone website.

It exposes ONLY read endpoints:
    GET /                         the viewer page
    GET /api/routes               list/search routes (compound / doi / author / year)
    GET /api/routes/{name}        one route: background + reactions + bipartite graph
    GET /api/mol.png              a molecule image (SMILES)
    GET /api/mol_mb.png           a tree-aligned molecule image (cached molblock)
    GET /api/reaction.png         a mapped-reaction image
    GET /api/whoami               login state (only if a login gate is configured)

There is NO editing, staging, review, hand-off, remapping, connectivity-check, or
raw-JSON endpoint here — data/routes is only ever read.

Run from the repo root:
    uvicorn scripts.dashboard.to_view.app:app --reload
"""

import base64
import collections
import functools
import hashlib
import io
import os
import re
import secrets
import threading
import time

from fastapi import FastAPI, HTTPException, Request
from fastapi.responses import HTMLResponse, Response
from fastapi.staticfiles import StaticFiles

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import rdDepictor

# Prefer CoordGen for cleaner 2D layouts (mirrors the edit build).
rdDepictor.SetPreferCoordGen(True)

from scripts.dashboard.to_view import chem

# --- paths -----------------------------------------------------------------
HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
# Original (curated) folder is configurable via DASH_ROUTES_DIR so the viewer can
# point at data/routes, data/cbyd_extra_batch, etc. Relative paths resolve
# against the repo root. (run.sh exposes a --routes flag.)
_routes_env = os.environ.get("DASH_ROUTES_DIR")
ROUTES_DIR = (
    os.path.join(REPO_ROOT, "data", "routes") if not _routes_env
    else _routes_env if os.path.isabs(_routes_env)
    else os.path.join(REPO_ROOT, _routes_env)
)
STATIC_DIR = os.path.join(HERE, "static")
print(f"[dashboard-view] original dir: {ROUTES_DIR}")


# --- helpers ---------------------------------------------------------------
def _safe_name(name):
    """Validate a route name (a bare filename, no extension); block path traversal."""
    if name.endswith(".json"):
        name = name[: -len(".json")]
    if (not name or "/" in name or "\\" in name or "\x00" in name
            or name in (".", "..")):
        raise HTTPException(status_code=400, detail=f"Invalid route name: {name!r}")
    base = os.path.realpath(ROUTES_DIR)
    full = os.path.realpath(os.path.join(ROUTES_DIR, f"{name}.json"))
    if os.path.commonpath([base, full]) != base:
        raise HTTPException(status_code=400, detail=f"Invalid route name: {name!r}")
    return name


def _routes_path(name):
    return os.path.join(ROUTES_DIR, f"{name}.json")


def _load_json(path):
    import json
    with open(path, "r") as f:
        return json.load(f)


def _load_route(name):
    name = _safe_name(name)
    path = _routes_path(name)
    if not os.path.exists(path):
        raise HTTPException(status_code=404, detail=f"Route not found: {name}")
    return _load_json(path)


def _route_numbers():
    """Map each route name -> its 1-based global alphabetical number (#1…#N)."""
    names = sorted((f[: -len(".json")] for f in os.listdir(ROUTES_DIR)
                    if f.endswith(".json")),
                   key=lambda n: (n.lower(), n))   # case-insensitive by first letter
    return {n: i + 1 for i, n in enumerate(names)}


# --- search index ----------------------------------------------------------
# background sits at the top of each file, so we read only the head and regex the
# fields rather than parse the JSON. `"name"` needs a quote right before `name`, so
# it does NOT also match `"drug name"` (captured separately).
_BG_NAME_RE = re.compile(r'"name"\s*:\s*"((?:[^"\\]|\\.)*)"')
_BG_DRUGNAME_RE = re.compile(r'"drug name"\s*:\s*"((?:[^"\\]|\\.)*)"')
_BG_AUTHOR_RE = re.compile(r'"author"\s*:\s*"((?:[^"\\]|\\.)*)"')
_BG_YEAR_RE = re.compile(r'"year"\s*:\s*"?([0-9]{3,4})')   # quoted or bare year
_DOI_RE = re.compile(r'"doi"\s*:\s*"((?:[^"\\]|\\.)*)"')

# The index caches by (dir mtime, file count) bucketed into 2s windows, so
# out-of-process edits show up within a couple of seconds without a restart.
_META_INDEX = {"key": None, "map": {}}


def _index_key():
    try:
        base = (os.path.getmtime(ROUTES_DIR), len(os.listdir(ROUTES_DIR)))
    except OSError:
        base = None
    return (base, int(time.time() // 2))


def _meta_index():
    """{name -> {compound, doi, author, year}} (all lowercased) for the sidebar
    search. One head-read per file covers all four fields."""
    key = _index_key()
    if _META_INDEX["key"] != key:
        m = {}
        for f in os.listdir(ROUTES_DIR):
            if not f.endswith(".json"):
                continue
            nm = f[: -len(".json")]
            try:
                with open(os.path.join(ROUTES_DIR, f)) as fh:
                    head = fh.read(2000)
                name = _BG_NAME_RE.search(head)
                drug = _BG_DRUGNAME_RE.search(head)
                doi = _DOI_RE.search(head)
                auth = _BG_AUTHOR_RE.search(head)
                year = _BG_YEAR_RE.search(head)
                compound = " | ".join(g.group(1) for g in (name, drug) if g).lower()
                m[nm] = {
                    "compound": compound,
                    "doi": (doi.group(1) if doi else "").lower(),
                    "author": (auth.group(1) if auth else "").lower(),
                    "year": (year.group(1) if year else "").lower(),
                }
            except Exception:
                m[nm] = {"compound": "", "doi": "", "author": "", "year": ""}
        _META_INDEX["key"] = key
        _META_INDEX["map"] = m
    return _META_INDEX["map"]


def _reactions_detail(synthesis):
    """Per-reaction detail (index-addressed) used by the route view."""
    out = []
    for i, rxn in enumerate(synthesis):
        out.append({
            "index": i,
            "reaction_id": str(rxn.get("reaction_id", i)),
            "number": i + 1,
            "raw_description": rxn.get("raw_description", "") or "",
            "reaction_smiles": rxn.get("reaction_smiles", "") or "",
            "full_reaction_smiles": rxn.get("full_reaction_smiles", "") or "",
            "full_reaction_smiles_mapped": rxn.get("full_reaction_smiles_mapped", "") or "",
            "mapping_confidence": rxn.get("mapping_confidence"),
            "parsed_reaction": rxn.get("parsed_reaction", {}) or {},
            "reagents": chem.reagent_summary(rxn),
        })
    return out


# --- images -----------------------------------------------------------------
def _png_response(png_bytes):
    return Response(content=png_bytes, media_type="image/png")


def _placeholder_png(text="invalid", size=(220, 180)):
    """A simple gray placeholder image for SMILES that RDKit cannot parse."""
    from PIL import Image, ImageDraw

    img = Image.new("RGB", size, (245, 245, 245))
    d = ImageDraw.Draw(img)
    d.rectangle([0, 0, size[0] - 1, size[1] - 1], outline=(200, 200, 200))
    d.text((size[0] // 2 - 24, size[1] // 2 - 6), text, fill=(150, 150, 150))
    buf = io.BytesIO()
    img.save(buf, format="PNG")
    return buf.getvalue()


@functools.lru_cache(maxsize=4096)
def _mol_png(smiles, size=220, color=False):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return _placeholder_png("invalid", (size, size))
    drawer = rdMolDraw2D.MolDraw2DCairo(size, size)
    opts = drawer.drawOptions()
    opts.padding = 0.08
    if not color:
        opts.useBWAtomPalette()      # monochrome: all atoms/bonds black
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    return drawer.GetDrawingText()


# chem.build_network attaches a per-node `molblock` (2D coords aligned to the
# reaction product). We stash those by content hash so the browser can fetch the
# already-oriented depiction with a short URL, then render from the stored coords.
MB_CACHE = collections.OrderedDict()
MB_CACHE_MAX = 20000


def _cache_graph_molblocks(graph):
    """Move each molecule node's `molblock` into MB_CACHE keyed by content hash,
    replacing it with a small `mb` key the client requests via /api/mol_mb.png."""
    for n in graph.get("nodes", []):
        mb = n.pop("molblock", None)
        if not mb:
            continue
        wedges = tuple(tuple(int(x) for x in p) for p in (n.pop("wedges", []) or []))
        key = hashlib.md5((mb + repr(wedges)).encode()).hexdigest()[:16]
        if key in MB_CACHE:
            MB_CACHE.move_to_end(key)
        else:
            MB_CACHE[key] = (mb, wedges)
            if len(MB_CACHE) > MB_CACHE_MAX:
                MB_CACHE.popitem(last=False)
        n["mb"] = key
    return graph


def _network(synthesis):
    """build_network + cache the per-node aligned molblocks (single call site)."""
    return _cache_graph_molblocks(chem.build_network(synthesis))


@functools.lru_cache(maxsize=8192)
def _molblock_png(mb, w=220, h=0, abbr=False, wedges=(), color=False):
    h = h or w
    mol = Chem.MolFromMolBlock(mb)                 # draw from the stored 2D coords
    if mol is None:
        return _placeholder_png("invalid", (w, h))
    if wedges:                                      # consistent wedge/dash across the tree
        chem.apply_wedges(mol, wedges)
    if abbr:                                        # compact depiction: Bn, Boc, PEG6, (CH2)n…
        amol, _ = chem.abbreviate_mol(mol)
        if amol is not None:
            mol = amol                              # keep the molblock's aligned coords
    drawer = rdMolDraw2D.MolDraw2DCairo(w, h)
    opts = drawer.drawOptions()
    opts.padding = 0.08
    if not color:
        opts.useBWAtomPalette()                     # monochrome: all atoms/bonds black
    if wedges:                                      # honour our wedging, don't recompute it
        mol = rdMolDraw2D.PrepareMolForDrawing(mol, wedgeBonds=False, addChiralHs=False)
        drawer.DrawMolecule(mol)
    else:
        rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol)
    drawer.FinishDrawing()
    return drawer.GetDrawingText()


# CoordGen (preferred globally) renders some REACTIONS as a blank canvas, so
# reaction drawing flips it off. Process-global flag, hence the lock.
_RXN_DRAW_LOCK = threading.Lock()


@functools.lru_cache(maxsize=2048)
def _reaction_png(rxn_smiles, sub_img=(380, 300), color=False):
    try:
        rxn = AllChem.ReactionFromSmarts(rxn_smiles, useSmiles=True)
        if rxn is None:
            raise ValueError("could not parse reaction")
        n_main = rxn.GetNumReactantTemplates() + rxn.GetNumProductTemplates()
        width = min(4000, max(2, n_main + 1) * sub_img[0])   # +1 leaves room for arrow/agents
        with _RXN_DRAW_LOCK:
            rdDepictor.SetPreferCoordGen(False)              # CoordGen blanks some reactions
            try:
                drawer = rdMolDraw2D.MolDraw2DCairo(width, sub_img[1])
                if not color:
                    drawer.drawOptions().useBWAtomPalette()      # monochrome
                drawer.DrawReaction(rxn)
                drawer.FinishDrawing()
                return drawer.GetDrawingText()
            finally:
                rdDepictor.SetPreferCoordGen(True)
    except Exception:
        return _placeholder_png("invalid reaction", (max(sub_img[0], 400), sub_img[1]))


# --- app -------------------------------------------------------------------
app = FastAPI(title="Reaction Route Viewer")

# Optional login gate (same scheme as the edit build), for a shared deployment.
#   DASH_USERS="alice:secret1,bob:secret2"  → named accounts
#   DASH_PASSCODE="shared"                   → one shared passcode
# If neither is set, the viewer is open.
def _parse_users(s):
    users = {}
    for pair in s.split(","):
        pair = pair.strip()
        if ":" in pair:
            u, p = pair.split(":", 1)
            if u.strip():
                users[u.strip()] = p
    return users


USERS = _parse_users(os.environ.get("DASH_USERS", ""))
PASSCODE = os.environ.get("DASH_PASSCODE", "").strip()


def _basic_creds(request):
    auth = request.headers.get("authorization", "")
    if not auth.startswith("Basic "):
        return None
    try:
        decoded = base64.b64decode(auth[6:]).decode("utf-8", "ignore")
    except Exception:
        return None
    return decoded.split(":", 1) if ":" in decoded else ["", decoded]


def _check_auth(request):
    """(ok, username) — username is set only when a named account matched."""
    creds = _basic_creds(request)
    if USERS:
        if not creds:
            return (False, None)
        u, p = creds
        exp = USERS.get(u)
        if exp is not None and secrets.compare_digest(p, exp):
            return (True, u)
        return (False, None)
    if PASSCODE:
        if not creds:
            return (False, None)
        return (secrets.compare_digest(creds[1], PASSCODE), None)
    return (True, None)   # open


@app.middleware("http")
async def _gate_and_no_cache(request, call_next):
    ok, user = _check_auth(request)
    if not ok:
        return Response("Authentication required", status_code=401,
                        headers={"WWW-Authenticate": 'Basic realm="Route viewer"'})
    request.state.auth_user = user
    response = await call_next(request)
    response.headers["Cache-Control"] = "no-store, must-revalidate"
    return response


@app.get("/api/whoami")
def whoami(request: Request):
    return {"user": getattr(request.state, "auth_user", None),
            "accounts": bool(USERS)}


@app.get("/", response_class=HTMLResponse)
def index():
    with open(os.path.join(STATIC_DIR, "index.html")) as f:
        return HTMLResponse(f.read())


@app.get("/api/routes")
def list_routes(q: str = "", num: str = "", bg: str = "", doi: str = "",
                author: str = "", year: str = "", limit: int = 0):
    """List routes with light metadata. Search by compound name (``bg``), DOI
    (``doi``), author (``author``), or year (``year``) — each an ANDed substring
    filter. ``q`` (filename) and ``num`` (global number/range) are also accepted."""
    if not os.path.isdir(ROUTES_DIR):
        raise HTTPException(status_code=500, detail=f"Routes dir missing: {ROUTES_DIR}")

    number_of = _route_numbers()                 # global #, alphabetical
    all_names = sorted(number_of, key=lambda n: number_of[n])

    ql = q.lower().strip()
    bgq = bg.lower().strip()
    doiq = doi.lower().strip()
    aq = author.lower().strip()
    yq = year.lower().strip()
    idx = _meta_index() if (bgq or doiq or aq or yq) else {}

    # parse the number filter: single "42", or a range "lo-hi" / "lo-" / "-hi"
    nq = num.strip()
    lo = hi = None
    num_active = bool(nq)
    if nq:
        try:
            if "-" in nq:
                a, b = nq.split("-", 1)
                lo = int(a) if a.strip() else None
                hi = int(b) if b.strip() else None
            else:
                lo = hi = int(nq)
        except ValueError:
            lo = hi = None          # unparseable → matches nothing

    def _num_ok(n):
        if not num_active:
            return True
        if lo is None and hi is None:
            return False
        if lo is not None and n < lo:
            return False
        if hi is not None and n > hi:
            return False
        return True

    items = []
    for name in all_names:
        meta = idx.get(name, {})
        if ql and ql not in name.lower():
            continue
        if bgq and bgq not in meta.get("compound", ""):
            continue
        if doiq and doiq not in meta.get("doi", ""):
            continue
        if aq and aq not in meta.get("author", ""):
            continue
        if yq and yq not in meta.get("year", ""):
            continue
        if not _num_ok(number_of[name]):
            continue
        items.append(name)
    total = len(items)
    if limit and limit > 0:
        items = items[:limit]

    out = []
    for name in items:
        meta = {"name": name, "number": number_of[name],
                "title": name, "author": None, "year": None, "n_steps": None}
        try:
            data = _load_json(_routes_path(name))
            bgd = data.get("background", {}) or {}
            meta["title"] = bgd.get("name") or name
            meta["author"] = bgd.get("author")
            meta["year"] = bgd.get("year")
            meta["n_steps"] = len(data.get("synthesis", []) or [])
        except Exception:
            pass
        out.append(meta)
    return {"total": total, "count": len(out), "n_all": len(all_names), "routes": out}


@app.get("/api/routes/{name}")
def get_route(name: str):
    """Full route detail: background, per-reaction info, and bipartite graph."""
    name = _safe_name(name)
    data = _load_route(name)
    synthesis = data.get("synthesis", []) or []
    reactions = _reactions_detail(synthesis)
    graph = _network(synthesis)
    n_targets = sum(1 for n in graph["nodes"]
                    if n.get("type") == "molecule" and n.get("role") == "target")
    return {
        "name": name,
        "background": data.get("background", {}) or {},
        "n_targets": n_targets,
        "reactions": reactions,
        "graph": graph,
    }


@app.get("/api/mol.png")
def mol_png(smiles: str, size: int = 220, color: int = 0):
    size = max(80, min(size, 500))
    return _png_response(_mol_png(smiles, size, bool(color)))


@app.get("/api/mol_mb.png")
def mol_mb_png(key: str, size: int = 220, w: int = 0, h: int = 0, abbr: int = 0, color: int = 0):
    """Render a tree-node molecule from its cached, tree-aligned 2D coordinates."""
    w = w or size
    h = h or w
    w = max(40, min(w, 3600))
    h = max(40, min(h, 3600))
    entry = MB_CACHE.get(key)
    if entry is None:
        return _png_response(_placeholder_png("reload", (w, h)))
    MB_CACHE.move_to_end(key)
    mb, wedges = entry
    return _png_response(_molblock_png(mb, w, h, bool(abbr), wedges, bool(color)))


@app.get("/api/reaction.png")
def reaction_png(name: str = "", index: int = -1, smiles: str = "", sub: int = 300, color: int = 0):
    """Render the mapped reaction. Pass an explicit ``smiles``, or ``name`` +
    ``index`` to render the currently stored mapped SMILES."""
    sub = max(160, min(sub, 900))
    sub_img = (int(sub * 1.25), sub)
    rxn_smiles = smiles
    if not rxn_smiles and name and index >= 0:
        name = _safe_name(name)
        synthesis = _load_route(name).get("synthesis", []) or []
        if 0 <= index < len(synthesis):
            rxn = synthesis[index]
            rxn_smiles = (rxn.get("full_reaction_smiles_mapped")
                          or rxn.get("full_reaction_smiles")
                          or rxn.get("reaction_smiles") or "")
    if not rxn_smiles:
        return _png_response(_placeholder_png("no reaction smiles", (700, 240)))
    return _png_response(_reaction_png(rxn_smiles, sub_img, bool(color)))


# Static assets (vendored JS/CSS). Mounted last so it does not shadow /api.
app.mount("/static", StaticFiles(directory=STATIC_DIR), name="static")


if __name__ == "__main__":
    import uvicorn
    uvicorn.run("scripts.dashboard.to_view.app:app", host="127.0.0.1", port=8001, reload=False)
