"""
build_static.py — bake the route viewer into ONE self-contained HTML file.

The live dashboard (app.py) renders every molecule / reaction image on the server
with RDKit. This script removes the server: it runs the same ``chem.build_network``
over every route once, ships the result inside the page, and lets the browser draw
the structures with RDKit.js (RDKit compiled to WebAssembly).

Each molecule node is stored as a CXSMILES string that carries its tree-aligned 2D
coordinates, its wedge/dash choices (``wU``/``wD``) and, for the abbreviated
depiction, its group labels (``$OMe$``). RDKit.js draws those exactly as the server
did. Each route is gzipped + base64'd separately and only unpacked when opened.

Usage (from the repo root, in an env with RDKit + NetworkX):
    python dashboard/build_static.py                       # -> dashboard.html
    python dashboard/build_static.py -o docs/index.html    # e.g. for GitHub Pages
    python dashboard/build_static.py --limit 20            # quick test build
    python dashboard/build_static.py --cdn                 # load RDKit.js from a CDN
                                                           #   (smaller file, needs internet)

The issue / PR / e-mail links can be set here (--repo, --email) or edited later
in the CONFIG block near the top of the generated HTML.
"""

import argparse
import base64
import gzip
import json
import os
import re
import subprocess
import sys
import time
import urllib.request
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(HERE, ".."))
STATIC_DIR = os.path.join(HERE, "static")
VENDOR_DIR = os.path.join(STATIC_DIR, "vendor")
sys.path.insert(0, REPO_ROOT)

from rdkit import Chem, RDLogger  # noqa: E402

from dashboard import chem  # noqa: E402

RDLogger.DisableLog("rdApp.*")

RDKIT_JS_VERSION = "2025.3.4-1.0.0"
RDKIT_JS_CDN = f"https://cdn.jsdelivr.net/npm/@rdkit/rdkit@{RDKIT_JS_VERSION}/dist/"

# parsed_reaction keys the viewer shows (names + the paired SMILES definitions).
PARSED_KEYS = [
    "GENERAL_REAGENTS", "CATALYSTS", "LIGANDS", "SOLVENTS", "MODIFIERS",
    "TEMPERATURE", "YIELD",
    "GENERAL_REAGENTS_DEFINITIONS", "CATALYSTS_DEFINITIONS",
    "LIGANDS_DEFINITIONS", "SOLVENTS_DEFINITIONS",
]

# CXSMILES extensions we need: 2D coords, atom labels (abbreviations), wedging.
_CX_FLAGS = (Chem.CXSmilesFields.CX_COORDS
             | Chem.CXSmilesFields.CX_ATOM_LABELS
             | Chem.CXSmilesFields.CX_BOND_CFG)


_CX_COORDS_RE = re.compile(r"\|\(([^)]*)\)")
_NUM_RE = re.compile(r"-?\d+\.\d+")


def _round_coord(m):
    s = f"{float(m.group(0)):.2f}".rstrip("0").rstrip(".")
    return "0" if s in ("-0", "") else s


def _cx(mol):
    """Mol -> CXSMILES carrying its current 2D coords, atom labels and wedges.
    Coords are cut to 2 decimals (sub-pixel at any size we draw) — they are most
    of the file size."""
    ps = Chem.SmilesWriteParams()
    cx = Chem.MolToCXSmiles(mol, ps, int(_CX_FLAGS))
    return _CX_COORDS_RE.sub(
        lambda m: "|(" + _NUM_RE.sub(_round_coord, m.group(1)) + ")", cx, count=1)


def _wedged(molblock, wedges):
    """Rebuild a node's mol from its aligned molblock with OUR wedge plan applied
    (same as app._molblock_png), so the CXSMILES carries those wedges."""
    mol = Chem.MolFromMolBlock(molblock)
    if mol is None:
        return None
    for b in mol.GetBonds():                 # drop the molblock's default wedging
        if b.HasProp("_MolFileBondCfg"):
            b.ClearProp("_MolFileBondCfg")
        if b.HasProp("_MolFileBondStereo"):
            b.ClearProp("_MolFileBondStereo")
    chem.apply_wedges(mol, wedges)
    return mol


def _mol_node(n):
    out = {"id": n["id"], "t": "m", "s": n["smiles"], "r": n["role"]}
    mb = n.get("molblock")
    if mb:
        wedges = [tuple(int(x) for x in p) for p in (n.get("wedges") or [])]
        mol = _wedged(mb, wedges) if wedges else Chem.MolFromMolBlock(mb)
        if mol is not None:
            out["cx"] = _cx(mol)
            if wedges:
                out["wd"] = 1                # draw with OUR wedges, don't recompute
            if n.get("ext_w_a") is not None:
                amol, _ = chem.abbreviate_mol(mol)
                if amol is not None:
                    out["cxa"] = _cx(amol)
                    out["ewa"], out["eha"] = n["ext_w_a"], n["ext_h_a"]
        out["ew"], out["eh"] = n.get("ext_w"), n.get("ext_h")
        out["na"] = n.get("n_atoms")
    return out


def export_route(path):
    """One route file -> (index row, gzipped+base64 payload) or (name, error)."""
    name = os.path.basename(path)[: -len(".json")]
    try:
        with open(path) as f:
            data = json.load(f)
        bg = data.get("background", {}) or {}
        synthesis = data.get("synthesis", []) or []
        graph = chem.build_network(synthesis)

        reactions = []
        for i, rxn in enumerate(synthesis):
            parsed = rxn.get("parsed_reaction", {}) or {}
            reactions.append({
                "rid": str(rxn.get("reaction_id", i)),
                "d": rxn.get("raw_description", "") or "",
                "m": rxn.get("full_reaction_smiles_mapped", "") or "",
                "p": {k: parsed[k] for k in PARSED_KEYS if parsed.get(k)},
            })

        nodes = []
        for n in graph["nodes"]:
            if n["type"] == "molecule":
                nodes.append(_mol_node(n))
            else:
                rn = {"id": n["id"], "t": "r", "i": n["index"]}
                if n.get("invalid"):
                    rn["inv"] = 1
                nodes.append(rn)

        payload = {
            "bg": bg,
            "rx": reactions,
            "nodes": nodes,
            "edges": [[e["source"], e["target"]] for e in graph["edges"]],
            "warn": graph.get("warnings", []),
        }
        blob = gzip.compress(json.dumps(payload, separators=(",", ":")).encode(), 9)
        compound = " | ".join(str(v) for v in (bg.get("name"), bg.get("drug name")) if v)
        row = {
            "name": name,
            "title": bg.get("name") or name,
            "author": bg.get("author"),
            "year": bg.get("year"),
            "n_steps": len(synthesis),
            "q_compound": compound.lower(),
            "q_doi": str(bg.get("doi") or "").lower(),
        }
        return row, base64.b64encode(blob).decode()
    except Exception as e:  # keep going; report at the end
        return name, f"{type(e).__name__}: {e}"


def _fetch(url, dest):
    if not os.path.exists(dest):
        print(f"  downloading {url}")
        os.makedirs(os.path.dirname(dest), exist_ok=True)
        urllib.request.urlretrieve(url, dest)
    with open(dest, "rb") as f:
        return f.read()


def _read(path):
    with open(path, encoding="utf-8") as f:
        return f.read()


def _inline_js(src):
    # a literal "</script" inside inlined JS would end the <script> element early
    return src.replace("</script", "<\\/script")


def _default_repo():
    """owner/name from the git remote, if it points at GitHub."""
    try:
        url = subprocess.run(["git", "-C", REPO_ROOT, "remote", "get-url", "origin"],
                             capture_output=True, text=True, check=True).stdout.strip()
    except Exception:
        return ""
    m = re.search(r"github\.com[:/]([^/]+/[^/]+?)(?:\.git)?$", url)
    return m.group(1) if m else ""


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--routes", default=os.path.join(REPO_ROOT, "data", "routes"),
                    help="folder of route JSON files (default: data/routes)")
    ap.add_argument("-o", "--out", default=os.path.join(REPO_ROOT, "dashboard.html"),
                    help="output HTML file (default: dashboard.html)")
    ap.add_argument("--repo", default=_default_repo(),
                    help="GitHub owner/name for the issue + PR links (default: git remote)")
    ap.add_argument("--branch", default="main", help="branch the PR links edit (default: main)")
    ap.add_argument("--email", default="", help="address for the 'Email us' link (default: none)")
    ap.add_argument("--data-path", default="data/routes",
                    help="path of the route files inside the GitHub repo (for PR links)")
    ap.add_argument("--cdn", action="store_true",
                    help="load RDKit.js from a CDN instead of inlining it (~9 MB smaller)")
    ap.add_argument("--limit", type=int, default=0, help="only the first N routes (testing)")
    ap.add_argument("-j", "--jobs", type=int, default=os.cpu_count() or 1)
    args = ap.parse_args()

    files = sorted((f for f in os.listdir(args.routes) if f.endswith(".json")),
                   key=lambda f: (f.lower(), f))          # same numbering as app.py
    if args.limit:
        files = files[: args.limit]
    paths = [os.path.join(args.routes, f) for f in files]
    print(f"Exporting {len(paths)} routes with {args.jobs} workers…")

    t0 = time.time()
    index, blobs, failed = [], {}, []
    with Pool(args.jobs) as pool:
        for k, (row, blob) in enumerate(pool.imap(export_route, paths, chunksize=4), 1):
            if isinstance(row, dict):
                row["number"] = k
                index.append(row)
                blobs[row["name"]] = blob
            else:
                failed.append((row, blob))
            if k % 250 == 0 or k == len(paths):
                print(f"  {k}/{len(paths)}  ({time.time() - t0:.0f}s)")
    for name, err in failed:
        print(f"  ! skipped {name}: {err}")

    if args.cdn:
        rdkit_js = f'<script src="{RDKIT_JS_CDN}RDKit_minimal.js"></script>'
        rdkit_wasm = ""
    else:
        js = _fetch(RDKIT_JS_CDN + "RDKit_minimal.js",
                    os.path.join(VENDOR_DIR, f"rdkit-{RDKIT_JS_VERSION}", "RDKit_minimal.js"))
        wasm = _fetch(RDKIT_JS_CDN + "RDKit_minimal.wasm",
                      os.path.join(VENDOR_DIR, f"rdkit-{RDKIT_JS_VERSION}", "RDKit_minimal.wasm"))
        rdkit_js = f"<script>{_inline_js(js.decode())}</script>"
        rdkit_wasm = base64.b64encode(wasm).decode()

    config = {
        "repo": args.repo, "branch": args.branch, "email": args.email,
        "dataPath": args.data_path.strip("/"),
        "rdkitWasmUrl": RDKIT_JS_CDN + "RDKit_minimal.wasm" if args.cdn else "",
        "built": time.strftime("%Y-%m-%d"),
    }
    vendor = "\n".join(
        f"<script>{_inline_js(_read(os.path.join(VENDOR_DIR, f)))}</script>"
        for f in ("dagre.min.js", "cytoscape.min.js", "cytoscape-dagre.js"))
    as_json = lambda o: json.dumps(o, separators=(",", ":")).replace("</", "<\\/")  # noqa: E731

    html = _read(os.path.join(STATIC_DIR, "standalone.html"))
    for key, val in {
        "/*__CONFIG__*/": json.dumps(config, indent=2),
        "/*__CSS__*/": _read(os.path.join(STATIC_DIR, "style.css")),
        "<!--__VENDOR__-->": vendor,
        "<!--__RDKIT_JS__-->": rdkit_js,
        "__RDKIT_WASM__": rdkit_wasm,
        "__INDEX__": as_json(index),
        "__ROUTES__": as_json(blobs),
        "/*__APP_JS__*/": _inline_js(_read(os.path.join(STATIC_DIR, "standalone.js"))),
    }.items():
        assert key in html, f"template placeholder missing: {key}"
        html = html.replace(key, val, 1)

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w", encoding="utf-8") as f:
        f.write(html)
    mb = os.path.getsize(args.out) / 1e6
    print(f"Wrote {args.out}  ({len(index)} routes, {mb:.1f} MB, {time.time() - t0:.0f}s)")


if __name__ == "__main__":
    main()
