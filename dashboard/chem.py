"""
chem.py — chemistry + reaction-network helpers for the dashboard.

This module is a *standalone adaptation* of the relevant pieces of
``scripts/route_utils.plot_rxn_network`` (and its helpers). It is deliberately
kept independent of ``scripts.route_utils`` / ``scripts.update_rxn_smiles`` so
the web backend does NOT pull in the heavy RXNMapper / torch import chain at
startup. Only RDKit + NetworkX are required here.

Two differences from the original ``plot_rxn_network``:
  * The graph is **bipartite**: molecule nodes AND explicit reaction nodes
    (reactant -> reaction -> product), so a reaction node is clickable in the
    UI and can carry the mapped SMILES + reagent metadata.
  * Nothing is drawn here; we emit a plain dict (nodes/edges) that the browser
    renders with Cytoscape.js. Image rendering lives in ``app.py``.

The molecule-cleaning behaviour (salt/complex removal, shared-species removal,
spectator halide stripping) mirrors the original so the set of molecule nodes
matches what ``plot_rxn_network`` would produce.
"""

import re
import threading
from collections import Counter, defaultdict

import numpy as np
from rdkit import Chem
from rdkit import RDLogger
from rdkit.Chem import AllChem, rdAbbreviations, rdDepictor, rdMolAlign, rdMolDescriptors
from rdkit.Geometry import Point2D

RDLogger.DisableLog("rdApp.*")

# CoordGen gives cleaner 2D layouts of complex (bridged / fused / macrocyclic)
# scaffolds. Process-global; also set in app.py, repeated here so build_network
# behaves the same when chem.py is imported on its own.
rdDepictor.SetPreferCoordGen(True)

# RDKit's default abbreviation set is small (37 simple groups + 8 linkers) and does
# NOT include benzyl, phenyl, or the common protecting groups — so O-CH2-Ph etc. never
# condensed. Add those ourselves. Order matters: CondenseMolAbbreviations matches in
# list order, so the more-specific acyl/silyl/benzyl groups must come BEFORE plain Ph
# (otherwise Ph would swallow the ring of a Bn/Bz and block the fuller abbreviation).
_CUSTOM_ABBREVS = rdAbbreviations.ParseAbbreviations("\n".join([
    "Boc   *C(=O)OC(C)(C)C            Boc",
    "Cbz   *C(=O)OCc1ccccc1           Cbz",
    "Fmoc  *C(=O)OCC1c2ccccc2-c2ccccc21  Fmoc",
    "TBDPS *[Si](c1ccccc1)(c1ccccc1)C(C)(C)C  TBDPS",
    "TBS   *[Si](C)(C)C(C)(C)C        TBS",
    "TIPS  *[Si](C(C)C)(C(C)C)C(C)C   TIPS",
    "TMS   *[Si](C)(C)C               TMS",
    "Tf    *S(=O)(=O)C(F)(F)F         Tf",
    "Ts    *S(=O)(=O)c1ccc(C)cc1      Ts",
    "Ms    *S(=O)(=O)C                Ms",
    # Benzoyl only when heteroatom-anchored (OBz ester / NHBz amide). A C-anchored
    # benzoyl is a phenyl KETONE (e.g. an acetophenone aldol product, -CC(=O)Ph) and
    # must stay drawn out, not collapse to "Bz" — same rule as Ac (OAc/NHAc, not a
    # methyl ketone).
    "OBz   *OC(=O)c1ccccc1            OBz",
    "NHBz  *NC(=O)c1ccccc1            NHBz",
    "Bn    *Cc1ccccc1                 Bn",
    "Ph    *c1ccccc1                  Ph",
]))

# Drop the generic standalone "Ac" abbreviation: it matches an acetyl on ANY atom,
# so a methyl KETONE (C-C(=O)-C) would render as "-Ac". Keep the heteroatom-anchored
# forms ("OAc", "NHAc") — acetyl is only condensed when it's an ester/amide.
_DEFAULT_ABBREVS = [a for a in rdAbbreviations.GetDefaultAbbreviations()
                    if a.label != "Ac"]

# Group abbreviations for compact depictions (Bn, Boc, Ts, OMe, nBu, nDec, …).
# We intentionally do NOT abbreviate bridging alkyl chains -(CH2)n-: condensing a
# plain methylene chain to "(CH2)6" hides structure and looks lopsided when only one
# of two similar chains gets caught. So drop the custom (CH2)n linkers AND RDKit's
# default alkyl chain linkers (Dec/Non/Oct/Hept), keeping only the PEG (oxygen)
# linkers. Terminal alkyl GROUPS (nBu, nDec, …) stay — those are in the defaults.
_PEG_LINKERS = [l for l in rdAbbreviations.GetDefaultLinkers()
                if l.label.startswith("PEG")]
_ABBREVS = list(_CUSTOM_ABBREVS) + _DEFAULT_ABBREVS + _PEG_LINKERS


def abbreviate_mol(mol):
    """Condense common groups / long chains into text labels. Returns
    (condensed_mol, sorted_labels); condensed_mol is None (labels []) if it fails
    or nothing condensed. The result carries `atomLabel` props (they do NOT survive
    a MolBlock roundtrip, so callers must render it directly)."""
    # A small molecule is already simple — abbreviating it (e.g. a whole reagent
    # collapsing to "OMe") only obscures it, so skip anything under 10 heavy atoms.
    if mol.GetNumHeavyAtoms() < 10:
        return None, []
    try:
        m = rdAbbreviations.CondenseMolAbbreviations(mol, _ABBREVS, maxCoverage=1.0)
    except Exception:
        return None, []
    if m is None or m.GetNumAtoms() >= mol.GetNumAtoms():
        return None, []
    # If the abbreviation ate the whole carbon skeleton (e.g. TBDPS-Cl -> "TBDPS·Cl",
    # TIPS-OTf -> "TIPS·OTf"), there's no backbone left to orient to — the preserved
    # coords can leave the label stacked vertically. Lay those out cleanly horizontal.
    if (not any(a.GetAtomicNum() == 6 for a in m.GetAtoms())
            and m.GetNumConformers()):
        _horizontalize(m)
    labels = sorted({a.GetProp("atomLabel") for a in m.GetAtoms()
                     if a.HasProp("atomLabel")})
    return m, labels


_TETRAHEDRAL = (Chem.ChiralType.CHI_TETRAHEDRAL_CW, Chem.ChiralType.CHI_TETRAHEDRAL_CCW)
# A shared fused ring system this big is a rigid core worth aligning to, even if a
# peripheral ring opens/closes in the reaction.
_FUSED_CORE_MIN = 10


def _wedge_plan(mol, inherited=()):
    """Which bond to draw as the wedge/dash for each stereocentre, as
    (stereocentre_idx, neighbour_idx) pairs. Starts from RDKit's default choice,
    then applies `inherited` overrides (pairs in THIS mol's atom indices) so a
    stereocentre shared with the parent is wedged on the SAME bond — giving a
    consistent depiction down the tree even when substituents change (e.g. the
    same centre showing a dash in the reactant and a wedge in the product). The
    up/down direction is computed at draw time (WedgeBond), so it always matches
    the actual configuration."""
    if mol.GetNumAtoms() == 0 or mol.GetNumConformers() == 0:
        return []
    tmp = Chem.Mol(mol)
    for b in tmp.GetBonds():
        b.SetBondDir(Chem.BondDir.NONE)
    Chem.WedgeMolBonds(tmp, tmp.GetConformer())     # default: begin atom = centre
    plan = {b.GetBeginAtomIdx(): b.GetEndAtomIdx()
            for b in tmp.GetBonds()
            if b.GetBondDir().name in ("BEGINWEDGE", "BEGINDASH")}
    # Only override centres that are ACTUAL tetrahedral stereocentres here: an atom
    # inherited from the parent may not be a stereocentre yet in this reactant (the
    # reaction creates it), and wedging a non-stereocentre crashes WedgeBond.
    chiral = {a.GetIdx() for a in mol.GetAtoms()
              if a.GetChiralTag() in _TETRAHEDRAL}
    for sc, nbr in inherited:                        # inherited bond wins
        if int(sc) in chiral and mol.GetBondBetweenAtoms(int(sc), int(nbr)) is not None:
            plan[int(sc)] = int(nbr)
    return [(sc, nbr) for sc, nbr in plan.items()]


def apply_wedges(mol, pairs):
    """Draw exactly `pairs` (stereocentre_idx, neighbour_idx) as wedge/dash on
    `mol` in place (clearing any existing wedging first). Direction is derived
    from the conformer + chirality by WedgeBond, so it's always correct. Callers
    then draw with wedgeBonds=False to honour these instead of re-wedging."""
    if not pairs or mol.GetNumConformers() == 0:
        return
    for b in mol.GetBonds():
        b.SetBondDir(Chem.BondDir.NONE)
    conf = mol.GetConformer()
    for sc, nbr in pairs:
        # Chirality can shift on the molblock roundtrip, so re-check here: WedgeBond
        # hard-asserts (crashing the whole render) if `sc` isn't a stereocentre.
        if mol.GetAtomWithIdx(int(sc)).GetChiralTag() not in _TETRAHEDRAL:
            continue
        bd = mol.GetBondBetweenAtoms(int(sc), int(nbr))
        if bd is not None:
            try:
                Chem.WedgeBond(bd, int(sc), conf)
            except Exception:
                pass


def _horizontalize(mol):
    """Rotate a molecule's 2D conformer so its LONGEST axis lies horizontal
    (principal-axis / PCA of the atom coordinates). Applied to the target root so
    the whole aligned tree inherits a horizontal orientation."""
    conf = mol.GetConformer()
    n = mol.GetNumAtoms()
    if n < 2:
        return
    pts = np.array([[conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y]
                    for i in range(n)])
    center = pts.mean(axis=0)
    p = pts - center
    # principal direction = eigenvector of the covariance with the largest eigenvalue
    evals, evecs = np.linalg.eigh(np.cov(p.T))
    vx, vy = evecs[:, int(np.argmax(evals))]
    theta = np.arctan2(vy, vx)                 # angle of the long axis
    c, s = np.cos(theta), np.sin(theta)
    R = np.array([[c, s], [-s, c]])            # rotate by -theta -> long axis to x
    rot = p @ R.T + center
    for i in range(n):
        conf.SetAtomPosition(i, (float(rot[i, 0]), float(rot[i, 1]), 0.0))


def _bbox(mol):
    """(width, height) of a molecule's 2D conformer bounding box."""
    conf = mol.GetConformer()
    xs = [conf.GetAtomPosition(i).x for i in range(mol.GetNumAtoms())]
    ys = [conf.GetAtomPosition(i).y for i in range(mol.GetNumAtoms())]
    if not xs:
        return 0.0, 0.0
    return max(xs) - min(xs), max(ys) - min(ys)


def _plain_coords(smiles):
    """A molecule with a plain CoordGen 2D layout, rotated so its longest axis is
    horizontal (used for the target root and as a fallback when there's nothing to
    align to). None if unparseable."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    rdDepictor.Compute2DCoords(mol)
    _horizontalize(mol)
    return mol


# Bridged/cage systems (≥ this many bridgehead atoms) are drawn with crossing
# bonds by flat CoordGen; for them we derive the 2D layout from a 3D conformer.
# But a large MACROCYCLE (ring ≥ this many atoms) is something CoordGen draws well
# and the 3D-derived layout crowds — so those stay on the CoordGen path even if they
# also contain a bridged core (e.g. eribulin: a 22-membered ring + polyether cage).
_BRIDGEHEAD_MIN = 3
_MACROCYCLE_MIN = 12
_BRIDGED_CACHE = {}          # canonical SMILES -> molblock (2D-from-3D), process-life


def _is_bridged(mol):
    try:
        if rdMolDescriptors.CalcNumBridgeheadAtoms(mol) < _BRIDGEHEAD_MIN:
            return False
        rings = mol.GetRingInfo().AtomRings()
        if rings and max(len(r) for r in rings) >= _MACROCYCLE_MIN:
            return False                            # big ring -> let CoordGen draw it
        return True
    except Exception:
        return False


def _bridged_coords(smiles):
    """A 2D layout derived from a 3D conformer, for bridged/cage molecules that flat
    CoordGen draws with crossing bonds (e.g. vinigrol, seychellene). Deterministic
    (seeded ETKDG) and cached by canonical SMILES, so the ~100 ms embed is paid once
    per unique molecule. Returns a mol with a 2D conformer, or None on failure."""
    canon = safe_canonicalize(smiles) or smiles
    if canon not in _BRIDGED_CACHE:
        # Some strained cages (e.g. seychellene) simply won't ETKDG-embed; cap the
        # attempt so a failure is quick, and cache BOTH success and failure so it's
        # paid at most once per molecule (a failure then falls back to flat CoordGen).
        result = None
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            try:
                mh = Chem.AddHs(mol)
                params = AllChem.ETKDGv3()
                params.randomSeed = 0xC0FFEE      # reproducible layout across rebuilds
                params.maxIterations = 50         # bound the (slow) failing embeds
                if AllChem.EmbedMolecule(mh, params) == 0:
                    AllChem.MMFFOptimizeMolecule(mh)
                    rdDepictor.GenerateDepictionMatching3DStructure(
                        mol, Chem.RemoveHs(mh), acceptFailure=True)
                    _horizontalize(mol)
                    result = Chem.MolToMolBlock(mol)
            except Exception:
                result = None
        _BRIDGED_CACHE[canon] = result or ""      # "" marks a failed / skipped embed
    mb = _BRIDGED_CACHE[canon]
    return Chem.MolFromMolBlock(mb) if mb else None


def _mapped_frags(mapped_smiles):
    """Parse an atom-mapped reaction SMILES into its reactant and product
    fragments, each keyed by canonical (map-stripped) SMILES -> the fragment mol
    WITH its atom-map numbers. Returns None when there's no usable mapping."""
    if not mapped_smiles:
        return None
    parts = split_rxn_smiles(mapped_smiles)
    if len(parts) < 3:
        return None

    def frags(block):
        out = {}
        for f in block.split("."):
            if not f:
                continue
            m = Chem.MolFromSmiles(f)
            if m is None or not any(a.GetAtomMapNum() for a in m.GetAtoms()):
                continue
            canon = safe_canonicalize(f, remove_map=True)
            if canon:
                out.setdefault(canon, m)
        return out

    return {"react": frags(parts[0]), "prod": frags(parts[-1])}


def _map_to_idx(mol, frag):
    """{atom_map_number: atom index in `mol`}, using `frag` (the same molecule,
    carrying map numbers) to look up which atom is which. None if they don't
    match structurally."""
    match = mol.GetSubstructMatch(frag)        # frag atom k -> mol atom index
    if not match:
        return None
    out = {}
    for k, mol_idx in enumerate(match):
        mn = frag.GetAtomWithIdx(k).GetAtomMapNum()
        if mn:
            out[mn] = mol_idx
    return out or None


# SetPreferCoordGen is a PROCESS-global flag; the pinning path below flips it off
# for a single layout, so serialize that window against other depiction.
_DEPICT_LOCK = threading.Lock()


def _layout_distorted(mol):
    """True if a 2D layout has a badly stretched/crammed bond OR a badly bent
    valence angle — the signature of pinning atoms onto coordinates that don't fit
    their geometry. A clean CoordGen layout has near-uniform bond lengths and
    near-ideal angles; good pins stay within a few degrees, bad ones bend ~20°+."""
    conf = mol.GetConformer()
    P = {i: np.array([conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y])
         for i in range(mol.GetNumAtoms())}
    bl = [np.linalg.norm(P[b.GetBeginAtomIdx()] - P[b.GetEndAtomIdx()])
          for b in mol.GetBonds()]
    if not bl:
        return False
    med = float(np.median(bl))
    if med <= 0:
        return True
    if (max(bl) / med > 1.6) or (med / min(bl) > 1.3):
        return True
    # A degree-2 atom bent far from its ideal angle = strained pin (e.g. a carbonyl
    # or ring CH2 squeezed onto the parent's coordinates).
    for a in mol.GetAtoms():
        if a.GetDegree() != 2:
            continue
        nb = [n.GetIdx() for n in a.GetNeighbors()]
        v1, v2 = P[nb[0]] - P[a.GetIdx()], P[nb[1]] - P[a.GetIdx()]
        n1, n2 = np.linalg.norm(v1), np.linalg.norm(v2)
        if n1 < 1e-6 or n2 < 1e-6:
            continue
        ang = np.degrees(np.arccos(np.clip(v1 @ v2 / (n1 * n2), -1, 1)))
        ideal = 120.0 if a.GetHybridization().name == "SP2" else 109.5
        if abs(ang - ideal) > 18:
            return True
    return False


def _clash_score(mol):
    """Count overlap problems in a 2D layout: non-bonded atoms sitting on top of
    each other, plus pairs of bonds that cross. Used to decide whether a layout
    needs a declash pass and whether that pass actually helped."""
    conf = mol.GetConformer()
    n = mol.GetNumAtoms()
    P = np.array([[conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y]
                  for i in range(n)])
    bonds = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in mol.GetBonds()]
    if not bonds:
        return 0
    L = float(np.median([np.linalg.norm(P[i] - P[j]) for i, j in bonds])) or 1.0
    bonded = {frozenset(b) for b in bonds}
    close = 0
    for i in range(n):
        for j in range(i + 1, n):
            if frozenset((i, j)) in bonded:
                continue
            if np.linalg.norm(P[i] - P[j]) < 0.6 * L:
                close += 1

    def _ccw(a, b, c):
        return (c[1] - a[1]) * (b[0] - a[0]) > (b[1] - a[1]) * (c[0] - a[0])

    cross = 0
    for x in range(len(bonds)):
        a, b = bonds[x]
        for y in range(x + 1, len(bonds)):
            c, d = bonds[y]
            if len({a, b, c, d}) < 4:                 # share an atom -> can't cross
                continue
            if (_ccw(P[a], P[c], P[d]) != _ccw(P[b], P[c], P[d]) and
                    _ccw(P[a], P[b], P[c]) != _ccw(P[a], P[b], P[d])):
                cross += 1
    return close + cross


def _relax_overlaps(mol, anchors, iters=60):
    """Minimally nudge a layout to remove overlaps: bond-length springs keep the
    scale, a soft repulsion pushes clashing non-bonded atoms apart, and a restraint
    pulls every atom back toward where it started — strongly for `anchors` (the
    atoms matched to the parent) so the alignment is preserved and only the free /
    flexible parts move. In place."""
    conf = mol.GetConformer()
    n = mol.GetNumAtoms()
    if n < 3:
        return
    P = np.array([[conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y]
                  for i in range(n)], float)
    P0 = P.copy()
    bonds = np.array([[b.GetBeginAtomIdx(), b.GetEndAtomIdx()]
                      for b in mol.GetBonds()])
    if len(bonds) == 0:
        return
    L = float(np.median(np.linalg.norm(P[bonds[:, 0]] - P[bonds[:, 1]], axis=1))) or 1.0
    bonded = np.zeros((n, n), bool)
    bonded[bonds[:, 0], bonds[:, 1]] = bonded[bonds[:, 1], bonds[:, 0]] = True
    np.fill_diagonal(bonded, True)
    anchor = np.full(n, 0.15)
    for i in anchors:
        anchor[i] = 3.0                               # matched atoms barely move
    R = 0.9 * L
    for _ in range(iters):
        F = (P0 - P) * anchor[:, None]                # restraint (keeps moves small)
        d = P[bonds[:, 0]] - P[bonds[:, 1]]           # bond springs
        dist = np.linalg.norm(d, axis=1)
        dist[dist < 1e-6] = 1e-6
        corr = ((dist - L) / dist * 0.5)[:, None] * d
        np.add.at(F, bonds[:, 0], -corr)
        np.add.at(F, bonds[:, 1], corr)
        diff = P[:, None, :] - P[None, :, :]          # pairwise repulsion for clashes
        dm = np.linalg.norm(diff, axis=2)
        mask = (dm < R) & (~bonded) & (dm > 1e-6)
        with np.errstate(invalid="ignore", divide="ignore"):
            mag = np.where(mask, (R - dm) / dm * 0.6, 0.0)
        F += np.nansum(mag[:, :, None] * diff, axis=1)
        P += np.clip(F, -0.3 * L, 0.3 * L) * 0.2
    for i in range(n):
        conf.SetAtomPosition(i, (float(P[i, 0]), float(P[i, 1]), 0.0))


def _declash(mol, anchors):
    """If a layout has overlaps, try a minimal relaxation (keeping `anchors` put)
    and accept it only if it reduces the overlap without distorting bonds — so it
    never makes a layout worse."""
    before = _clash_score(mol)
    if before == 0:
        return mol
    relaxed = Chem.Mol(mol)
    _relax_overlaps(relaxed, anchors)
    if _clash_score(relaxed) < before and not _layout_distorted(relaxed):
        return relaxed
    return mol


def _flip_to_match_parent(mol, pairs, parent_mol):
    """After a molecule is horizontalized, its long axis is horizontal but the
    PCA that chose that axis has an arbitrary sign, so the molecule may point
    left-or-right (and up-or-down) at random — making consecutive tree molecules
    flip-flop even when they share a scaffold. Pick, among the two chirality-safe
    orientations that keep the axis horizontal (identity and 180° rotation), the
    one whose shared atoms best match the parent's coordinates. 180° is a PROPER
    rotation (negate x AND y), so it preserves the 2D drawing's handedness — a
    mirror would flip wedge/dash. Mutates `mol`'s conformer in place."""
    if len(pairs) < 2:
        return
    conf = mol.GetConformer()
    pconf = parent_mol.GetConformer()
    # centre the parent's shared atoms so the comparison is translation-invariant
    # (horizontalize already centred `mol` on its own centroid).
    pcen_x = sum(pconf.GetAtomPosition(p).x for _, p in pairs) / len(pairs)
    pcen_y = sum(pconf.GetAtomPosition(p).y for _, p in pairs) / len(pairs)
    ccen_x = sum(conf.GetAtomPosition(c).x for c, _ in pairs) / len(pairs)
    ccen_y = sum(conf.GetAtomPosition(c).y for c, _ in pairs) / len(pairs)

    def sse(flip):
        s = 0.0
        for c, p in pairs:
            cp = conf.GetAtomPosition(c)
            dx = flip * (cp.x - ccen_x) - (pconf.GetAtomPosition(p).x - pcen_x)
            dy = flip * (cp.y - ccen_y) - (pconf.GetAtomPosition(p).y - pcen_y)
            s += dx * dx + dy * dy
        return s

    if sse(-1.0) < sse(1.0):                       # 180° rotation fits the parent better
        for i in range(mol.GetNumAtoms()):
            q = conf.GetAtomPosition(i)
            conf.SetAtomPosition(i, (2 * ccen_x - q.x, 2 * ccen_y - q.y, 0.0))


# Below this many heavy atoms a portrait layout is fine (a small fragment), so the
# landscape rescue only kicks in for sizable molecules.
_LANDSCAPE_MIN_ATOMS = 12


def _as_landscape(mol, pairs, parent_mol):
    """Aligning/pinning to a parent can hand back a cramped PORTRAIT layout (taller
    than wide) — and because each reactant is oriented to its product, one such
    layout cascades: every descendant pins to it and comes out narrower still. For
    anything but a small fragment a portrait depiction is essentially always an
    artifact (molecules read best laid out along a horizontal axis). Rotate it
    landscape (longest axis horizontal) and flip to match the parent so left/right
    stays consistent. Mutates and returns `mol`."""
    if mol.GetNumHeavyAtoms() >= _LANDSCAPE_MIN_ATOMS:
        w, h = _bbox(mol)
        if h > w * 1.2:
            _horizontalize(mol)
            _flip_to_match_parent(mol, pairs, parent_mol)
    return mol


def _n_shared_rings(mol, shared_idxs):
    """How many of `mol`'s rings are composed ENTIRELY of shared atoms. Comparing
    this between reactant and product tells whether a ring formed/broke among the
    shared scaffold (a real cyclization) vs a pendant ring added on new atoms."""
    return sum(1 for r in mol.GetRingInfo().AtomRings()
               if shared_idxs.issuperset(r))


def _largest_fused_shared_ring(mol, shared_idxs):
    """Size of the largest connected (fused) ring system built only from shared ring
    atoms. A big value means there's a large rigid core worth aligning to (align it
    even if a peripheral ring opens); a small value means the shared rings are small
    and isolated (a cyclization there is better shown as a stretched chain)."""
    ring_atoms = [a.GetIdx() for a in mol.GetAtoms()
                  if a.IsInRing() and a.GetIdx() in shared_idxs]
    adj = {a: set() for a in ring_atoms}
    for b in mol.GetBonds():
        if b.IsInRing():
            i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
            if i in adj and j in adj:
                adj[i].add(j)
                adj[j].add(i)
    seen, best = set(), 0
    for start in ring_atoms:
        if start in seen:
            continue
        stack, size = [start], 0
        while stack:
            x = stack.pop()
            if x in seen:
                continue
            seen.add(x)
            size += 1
            stack.extend(y for y in adj[x] if y not in seen)
        best = max(best, size)
    return best


def _aligned_child(smiles, react_frag, parent_map2idx, parent_mol):
    """Depict `smiles` (a reactant) oriented to the parent (product) via reaction
    atom-map numbers, so a shared scaffold keeps one orientation down the tree.

    Two strategies, preferring the least-distorting one:

    1. RIGID orient: give the molecule its own clean CoordGen layout and rotate /
       translate it (rdMolAlign) so the shared atoms best-fit the parent. This
       keeps a natural, uncontorted shape and works for most molecules (a rigid
       fused-ring scaffold is drawn the same way in both, so it fits exactly).

    2. PIN: only if (1) fits poorly — i.e. the shared part has rotatable bonds and
       CoordGen drew it as a different rotamer in each molecule, so no rigid
       rotation reconciles it (e.g. trityl / pyrrolidine). Then we pin the shared
       atoms onto the parent's coordinates with the classic depictor (CoordGen
       ignores a coordinate template) and lay the rest out around them.

    Preferring (1) avoids the pinned layout occasionally contorting or flipping a
    simple molecule (which makes its stereochemistry hard to read). Returns None
    only when the SMILES can't be parsed."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    rdDepictor.Compute2DCoords(mol)                  # its own natural layout
    child_idx = _map_to_idx(mol, react_frag)
    if not child_idx:
        return mol
    pairs = [(cidx, parent_map2idx[mn]) for mn, cidx in child_idx.items()
             if mn in parent_map2idx]
    if len(pairs) < 2:                               # need >=2 anchors for orientation
        return mol
    # 1) Rigidly orient the molecule's own clean layout onto the parent. This is the
    #    fallback; it never distorts, but two independent CoordGen layouts can differ
    #    (a shared ring drawn mirrored between reactant and product), so a rigid fit
    #    alone can leave the ring flipped.
    try:
        rdMolAlign.AlignMol(mol, parent_mol, atomMap=pairs, maxIters=50)
    except Exception:
        pass
    natural = Chem.Mol(mol)
    # Did a ring form / break AMONG THE SHARED ATOMS — a real cyclization or
    # ring-opening (e.g. a chain closes to a ring, or a monocycle becomes bicyclic)?
    # Count rings whose atoms are ALL shared, on each side: a difference means the
    # shared scaffold changes ring count. This is distinct from a pendant ring merely
    # ADDED on new atoms (e.g. a tosyl), which leaves the shared-ring count unchanged.
    child_shared = {c for c, _ in pairs}
    parent_shared = {p for _, p in pairs}
    ring_change = (_n_shared_rings(mol, child_shared)
                   != _n_shared_rings(parent_mol, parent_shared))
    # ...but if a LARGE fused ring core is shared, keep aligning to it (a peripheral
    # ring opening, e.g. a lactone off a big polycyclic, shouldn't unalign the core).
    big_core = _largest_fused_shared_ring(mol, child_shared) >= _FUSED_CORE_MIN
    ring_pairs = [(c, p) for c, p in pairs if mol.GetAtomWithIdx(c).IsInRing()]
    if (ring_change and not big_core) or len(ring_pairs) < 3:
        # A cyclization/ring-opening among small, isolated shared rings — aligning to
        # the product would coil the chain or rotate the molecule vertical. Lay it out
        # as a clean chain instead — longest axis horizontal. Then flip it (180°, a
        # chirality-safe proper rotation) if that better matches the parent, so the
        # stretched chain still points the same way as the molecule above it in the
        # tree instead of flip-flopping on the PCA sign.
        _horizontalize(natural)
        _flip_to_match_parent(natural, pairs, parent_mol)
        return natural
    # 2) Pin the shared atoms exactly onto the parent so a shared scaffold is drawn
    #    identically (no flip / drift). The guards below fall back to `natural` when
    #    pinning would distort or overlap the molecule.
    #    Anchor the ring core PLUS any shared stereocentre and its shared neighbours:
    #    an acyclic stereocentre hanging off a pinned ring (e.g. the benzylic carbon of
    #    a pyridyl amine) would otherwise be laid out freely and can come out mirrored
    #    between reactant and product, flipping its wedge/dash.
    pair_of = dict(pairs)
    anchor = {c for c, _ in ring_pairs}
    for a in mol.GetAtoms():
        if a.GetChiralTag() in _TETRAHEDRAL and a.GetIdx() in pair_of:
            anchor.add(a.GetIdx())
            anchor.update(n.GetIdx() for n in a.GetNeighbors()
                          if n.GetIdx() in pair_of)
    anchor_pairs = [(c, pair_of[c]) for c in anchor]
    pconf = parent_mol.GetConformer()
    coord_map = {cidx: Point2D(pconf.GetAtomPosition(pidx).x,
                               pconf.GetAtomPosition(pidx).y)
                 for cidx, pidx in anchor_pairs}
    with _DEPICT_LOCK:
        try:
            rdDepictor.SetPreferCoordGen(False)      # classic honors the template
            rdDepictor.Compute2DCoords(mol, coordMap=coord_map)
        except Exception:
            rdDepictor.Compute2DCoords(mol)
        finally:
            rdDepictor.SetPreferCoordGen(True)
    if _layout_distorted(mol):                       # pinning made it weird -> relax
        return _as_landscape(natural, pairs, parent_mol)  # clean shape, approx orient
    # Nudge out any overlap while keeping the anchored (aligned) atoms put.
    relaxed = _declash(mol, [cidx for cidx, _ in anchor_pairs])
    # If pinning leaves the layout NOTICEABLY more overlapping (e.g. an ester stacked
    # on a pyridine) than the clean natural layout, prefer natural — a clean,
    # approximately-oriented depiction beats an aligned-but-overlapping one. Use a
    # margin so a single inherent crossing (e.g. a crowded trityl) still stays pinned.
    if _clash_score(relaxed) >= _clash_score(natural) + 2:
        return _as_landscape(natural, pairs, parent_mol)
    # A pinned layout that came out cramped-portrait for a sizable molecule is an
    # orientation artifact that cascades down the branch — lay the clean CoordGen
    # shape out landscape instead (this re-orients the whole downstream chain, since
    # each reactant pins to its now-landscape product).
    rw, rh = _bbox(relaxed)
    if mol.GetNumHeavyAtoms() >= _LANDSCAPE_MIN_ATOMS and rh > rw * 1.2:
        return _as_landscape(natural, pairs, parent_mol)
    return relaxed

# Neutralizer for protonation-invariant matching (an acid and its conjugate base —
# HCl/[Cl-], TFA/TFA- — map to the same neutral form). Optional: falls back to
# exact matching if the standardizer isn't available.
try:
    from rdkit.Chem.MolStandardize import rdMolStandardize
    _UNCHARGER = rdMolStandardize.Uncharger()
except Exception:  # pragma: no cover - depends on RDKit build
    _UNCHARGER = None


def _protonation_key(smiles):
    """Canonical SMILES of the neutralized (uncharged) molecule, so a species and
    its protonation/charge variant share a key (e.g. 'Cl' and '[Cl-]' -> 'Cl';
    TFA and TFA- -> 'O=C(O)C(F)(F)F'). None if unparseable."""
    if not smiles:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    for a in mol.GetAtoms():
        a.SetAtomMapNum(0)
    if _UNCHARGER is not None:
        try:
            mol = _UNCHARGER.uncharge(mol)
        except Exception:
            pass
    try:
        return Chem.MolToSmiles(mol)
    except Exception:
        return None


# ---------------------------------------------------------------------------
# SMILES utilities (adapted from route_utils / update_rxn_smiles)
# ---------------------------------------------------------------------------
def split_rxn_smiles(rxn_smiles):
    """Split a reaction SMILES on '>' that is not part of '->' or '<-'."""
    return re.split(r"(?<!-)>(?!-)", rxn_smiles)


def canonicalize_smiles(smiles, remove_map=True):
    """Canonicalize a SMILES string, optionally stripping atom-map numbers.

    Raises ValueError on an unparseable SMILES (same contract as the original
    in update_rxn_smiles.canonicalize_smiles).
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")
    if remove_map:
        for atom in mol.GetAtoms():
            atom.SetAtomMapNum(0)
    try:
        smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
        mol = Chem.MolFromSmiles(smiles)
    except Exception:
        pass
    return Chem.MolToSmiles(mol, isomericSmiles=True)


def safe_canonicalize(smiles, remove_map=True):
    """canonicalize_smiles that returns None instead of raising."""
    try:
        return canonicalize_smiles(smiles, remove_map=remove_map)
    except Exception:
        return None


# Canonicalized forms of spectator salt-forming acids and their conjugate anions,
# used to detect/strip counter-ions that should not become network nodes/roots.
_HALOGEN_ACID_SMILES_RAW = [
    "Cl", "Cl-", "[Cl-]", "Cl[H]", "ClH", "Cl.H", "HCl", "[H][Cl]",
    "Br", "Br-", "[Br-]", "Br[H]", "BrH", "Br.H", "HBr", "[H][Br]",
    "I", "I-", "[I-]", "IH", "I.H", "HI", "[H][I]",
    "F", "F-", "[F-]", "FH", "F.H", "HF", "[H][F]",
]

# Common non-halide counter-ions (the acid and every anionic protonation state) —
# treated like halide counter-ions for spectator-salt detection/stripping.
_OXOACID_SALT_SMILES_RAW = [
    "OP(=O)(O)O", "OP(=O)(O)[O-]", "OP(=O)([O-])[O-]", "[O-]P(=O)([O-])[O-]",  # H3PO4 / H2PO4- / HPO4-2 / PO4-3
    "OS(=O)(=O)O", "OS(=O)(=O)[O-]", "[O-]S(=O)(=O)[O-]",                       # H2SO4 / HSO4- / SO4-2
    "Cc1ccc(S(=O)(=O)O)cc1", "Cc1ccc(S(=O)(=O)[O-])cc1",                        # TsOH / TsO- (tosylate)
    "OS(C)(=O)=O", "[O-]S(C)(=O)=O",                                            # MsOH / MsO- (mesylate)
    "OS(=O)(=O)C(F)(F)F", "[O-]S(=O)(=O)C(F)(F)F",                              # TfOH / TfO- (triflate)
    "OC(=O)C(F)(F)F", "[O-]C(=O)C(F)(F)F",                                      # TFA / TFA- (trifluoroacetate)
]


def get_halogen_acid_set():
    """Canonical SMILES of spectator counter-ions (halide acids/anions plus the
    phosphate & sulfate oxo-acids and their anions)."""
    acids = set()
    for salt in _HALOGEN_ACID_SMILES_RAW + _OXOACID_SALT_SMILES_RAW:
        can = safe_canonicalize(salt)
        if can:
            acids.add(can)
    return acids


def remove_halide_salts(rxn_smiles):
    """Strip spectator halide species (HI/[I-]/HCl/...) from a side only when a
    real organic molecule (>1 carbon) still remains on that side."""
    halogen_acids = get_halogen_acid_set()

    def count_carbons(smi):
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            return 0
        return sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 6)

    def strip_side(side_smi):
        comps = [c for c in side_smi.split(".") if c]
        non_halide, halide = [], []
        for c in comps:
            can = safe_canonicalize(c)
            if can is not None and can in halogen_acids:
                halide.append(c)
            else:
                non_halide.append(c)
        if halide and any(count_carbons(c) > 1 for c in non_halide):
            return non_halide
        return comps

    parts = split_rxn_smiles(rxn_smiles)
    if len(parts) < 3:
        return rxn_smiles
    reactants_smi, agents, products_smi = parts[0], parts[1], parts[-1]
    reactants_new = strip_side(reactants_smi)
    products_new = strip_side(products_smi)
    return ".".join(reactants_new) + ">" + agents + ">" + ".".join(products_new)


def _comp_map_nums(comp):
    """Atom-map numbers used in a single SMILES component (empty if unparseable)."""
    mol = Chem.MolFromSmiles(comp)
    if mol is None:
        return set()
    return {a.GetAtomMapNum() for a in mol.GetAtoms() if a.GetAtomMapNum()}


def remove_spectator_halides(rxn_smiles):
    """Mapping-aware halide cleanup for the tree.

    A halide (HBr/HCl/[Cl-]/…) is KEPT only if it is *incorporated into an organic
    product* — i.e. one of its atom-map numbers appears in a product component that
    is itself not a halide salt (e.g. HBr → [Br:1][CH:2]… keeps HBr). A halide that
    maps only to another halide/counter-ion (HCl → [Cl-]) or isn't mapped in is a
    spectator and is stripped (when a >1-carbon organic remains on the side). With
    no atom maps at all, nothing is "incorporated", so this matches the plain
    remove_halide_salts behaviour."""
    halogen_acids = get_halogen_acid_set()

    parts = split_rxn_smiles(rxn_smiles)
    if len(parts) < 3:
        return rxn_smiles
    reactants_smi, agents, products_smi = parts[0], parts[1], parts[-1]

    # atom-map numbers that belong to a NON-halide (organic) product component
    organic_map_nums = set()
    for c in products_smi.split("."):
        if not c:
            continue
        can = safe_canonicalize(c)
        if can is not None and can in halogen_acids:
            continue                       # a standalone halide product doesn't count
        organic_map_nums |= _comp_map_nums(c)

    def count_carbons(smi):
        mol = Chem.MolFromSmiles(smi)
        return 0 if mol is None else sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 6)

    def strip_side(side_smi):
        comps = [c for c in side_smi.split(".") if c]
        non_halide, halide = [], []
        for c in comps:
            can = safe_canonicalize(c)
            is_halide = can is not None and can in halogen_acids
            incorporated = bool(_comp_map_nums(c) & organic_map_nums)
            if is_halide and not incorporated:
                halide.append(c)           # spectator halide -> candidate to strip
            else:
                non_halide.append(c)       # organic, or a halide fixed into the product
        if halide and any(count_carbons(c) > 1 for c in non_halide):
            return non_halide
        return comps

    reactants_new = strip_side(reactants_smi)
    products_new = strip_side(products_smi)
    return ".".join(reactants_new) + ">" + agents + ">" + ".".join(products_new)


def check_if_reactant_product_same(rxn_smiles):
    """True only if reactants == products (a genuine salt/complex no-op like
    ``A.HCl >> A.HCl``), ignoring spectator counter-ions (halide/phosphate/sulfate/
    sulfonate/TFA) when both sides carry them alongside other species. Comparison is
    charge-SENSITIVE, so a real protonation step (``A >> AH+``) is NOT dropped — it
    is shown as its own node (e.g. R-NO2 -> R-NH2 -> R-NH3+)."""
    salts = get_halogen_acid_set()

    def split_salt(smi):
        comps = [c for c in smi.split(".") if c]
        canon = [(c, safe_canonicalize(c)) for c in comps]
        return ([c for c, ca in canon if ca not in salts],
                [c for c, ca in canon if ca in salts])

    parts = split_rxn_smiles(rxn_smiles)
    if len(parts) < 3:
        return False
    reactant_smiles, product_smiles = parts[0], parts[2]

    nh_react, h_react = split_salt(reactant_smiles)
    nh_prod, h_prod = split_salt(product_smiles)

    if h_react and h_prod and nh_react and nh_prod:
        rc, pc = ".".join(nh_react), ".".join(nh_prod)
    else:
        rc, pc = reactant_smiles, product_smiles
    a, b = safe_canonicalize(rc), safe_canonicalize(pc)
    return a is not None and a == b


def clean_rxn_duplicates(rxn_smiles):
    """Remove molecules unchanged between both sides (e.g. A.B>>C.B -> A>>C)."""
    parts = split_rxn_smiles(rxn_smiles)
    if len(parts) < 3:
        return rxn_smiles
    reactants_smi, agents, products_smi = parts[0], parts[1], parts[-1]
    reactants = reactants_smi.split(".")
    products = products_smi.split(".")

    reactants_canon = [safe_canonicalize(r) for r in reactants]
    products_canon = [safe_canonicalize(p) for p in products]
    common = set(reactants_canon) & set(products_canon)
    common.discard(None)

    reactants_new = [
        r for r, rc in zip(reactants, reactants_canon)
        if rc not in common or reactants_canon.count(rc) > products_canon.count(rc)
    ]
    products_new = [
        p for p, pc in zip(products, products_canon)
        if pc not in common or products_canon.count(pc) > reactants_canon.count(pc)
    ]

    # Also drop a spectator that appears on BOTH sides differing only by protonation
    # / charge (e.g. HCl / [Cl-], TFA / TFA-). Guarded: only when both sides keep >1
    # component, so a genuine sole reactant/product is never stripped.
    if len(reactants_new) > 1 and len(products_new) > 1:
        r_keys = [_protonation_key(r) for r in reactants_new]
        p_keys = [_protonation_key(p) for p in products_new]
        shared = {k for k in r_keys if k} & {k for k in p_keys if k}
        if shared:
            r2 = [r for r, k in zip(reactants_new, r_keys) if k not in shared]
            p2 = [p for p, k in zip(products_new, p_keys) if k not in shared]
            if r2 and p2:                 # never empty a side
                reactants_new, products_new = r2, p2

    return ".".join(reactants_new) + ">" + agents + ">" + ".".join(products_new)


# ---------------------------------------------------------------------------
# Reagent summary
# ---------------------------------------------------------------------------
_REAGENT_CATEGORIES = [
    ("GENERAL_REAGENTS", "Reagents"),
    ("CATALYSTS", "Catalysts"),
    ("LIGANDS", "Ligands"),
    ("SOLVENTS", "Solvents"),
    ("MODIFIERS", "Modifiers"),
    ("TEMPERATURE", "Temperature"),
    ("YIELD", "Yield"),
]


def reagent_summary(reaction):
    """Build a structured, display-friendly reagent summary for a reaction
    dict from a route JSON's ``synthesis`` list."""
    parsed = reaction.get("parsed_reaction", {}) or {}
    out = []
    for key, label in _REAGENT_CATEGORIES:
        values = parsed.get(key) or []
        if values:
            out.append({"label": label, "key": key, "values": list(values)})
    return out


# ---------------------------------------------------------------------------
# Bipartite reaction-network builder
# ---------------------------------------------------------------------------
def _clean_reaction_for_nodes(rxn_smiles):
    """Apply the same cleaning chain as plot_rxn_network and return
    (reactant_canon_list, product_canon_list) of *molecule node* SMILES, or
    None if the reaction should be skipped (salt/complex formation) or invalid."""
    if check_if_reactant_product_same(rxn_smiles):
        return None
    cleaned = clean_rxn_duplicates(rxn_smiles)
    # mapping-aware: keep a halide that's incorporated into an organic product
    # (e.g. HBr adding across a diene), strip spectator counter-ions (HCl -> [Cl-]).
    cleaned = remove_spectator_halides(cleaned)
    parts = split_rxn_smiles(cleaned)
    if len(parts) < 3:
        return None
    reactants = [safe_canonicalize(r) for r in parts[0].split(".") if r]
    products = [safe_canonicalize(p) for p in parts[-1].split(".") if p]
    reactants = [r for r in reactants if r]
    products = [p for p in products if p]
    if not reactants or not products:
        return None
    return reactants, products


def _reaction_cycles(per_rxn):
    """Find cycles in the reaction-flow graph (a product of one reaction is a
    reactant of another). A cycle means a molecule is regenerated from something
    downstream of it, so nothing is left as a terminal product — the usual cause
    of "0 targets". Returns short human strings like "Rxn 3 ⇄ Rxn 4 (…)".

    ``per_rxn[i]`` is ``(reactants, products)`` (canonical SMILES) or ``None``.
    """
    n = len(per_rxn)
    prods = [set(c[1]) if c else set() for c in per_rxn]
    reacts = [set(c[0]) if c else set() for c in per_rxn]
    adj = {i: {j for j in range(n) if j != i and prods[i] & reacts[j]}
           for i in range(n)}

    # Tarjan's strongly-connected components; any SCC with >1 reaction is a cycle.
    idx = {}; low = {}; on = set(); stack = []; sccs = []; c = [0]
    def strong(v):
        idx[v] = low[v] = c[0]; c[0] += 1
        stack.append(v); on.add(v)
        for w in adj[v]:
            if w not in idx:
                strong(w); low[v] = min(low[v], low[w])
            elif w in on:
                low[v] = min(low[v], idx[w])
        if low[v] == idx[v]:
            comp = []
            while True:
                w = stack.pop(); on.discard(w); comp.append(w)
                if w == v:
                    break
            sccs.append(comp)
    import sys as _sys
    _lim = _sys.getrecursionlimit()
    _sys.setrecursionlimit(max(_lim, n * 4 + 100))
    try:
        for v in range(n):
            if v not in idx:
                strong(v)
    finally:
        _sys.setrecursionlimit(_lim)

    out = []
    for comp in sccs:
        if len(comp) < 2:
            continue
        members = sorted(comp)
        mset = set(members)
        # a molecule made inside the loop and consumed inside the loop
        shared = ""
        for i in members:
            common = prods[i] & set().union(*(reacts[j] for j in mset if j != i))
            if common:
                shared = sorted(common, key=len)[0]
                break
        nums = [i + 1 for i in members]
        join = " ⇄ " if len(nums) == 2 else " → "
        label = join.join(f"Rxn {x}" for x in nums)
        if len(nums) > 2:
            label += f" → Rxn {nums[0]}"
        if shared:
            s = shared if len(shared) <= 44 else shared[:41] + "…"
            label += f"  ({s} is remade from its own product)"
        out.append(label)
    return out


def build_network(synthesis, align=True):
    """Build a bipartite reaction network from a route's ``synthesis`` list.

    ``align=False`` skips the (expensive) 2D depiction / atom-map alignment work,
    producing the same topology (roles, edges, warnings) with no ``molblock`` on
    the nodes — used by the batch connectivity scan.

    Returns a dict:
        {
          "nodes": [
             {"id","type":"molecule","smiles","role"} |
             {"id","type":"reaction","index","reaction_id","number","label",
              "raw_description","reagents":[...],
              "has_mapped":bool, "mapping_confidence":float|None}
          ],
          "edges": [{"source","target"}],
          "warnings": [str, ...],
          "n_steps": int,
        }

    Topology — one molecule node per *occurrence* (no molecule points to more
    than one reaction):
      * One reaction node per synthesis entry (keyed by ``index`` — reaction_id
        is NOT unique within a route; convergent syntheses reuse 1,2,3… per
        branch).
      * Every reactant of a reaction is its OWN molecule node. If that molecule
        was made by an earlier reaction, we connect ``producer -> molecule ->
        this reaction`` — so a building block used in N steps appears as N copies
        (N branches from its producer), instead of one shared node fanning out to
        N reactions. A reactant that nothing produces is a starting material.
      * A product that is consumed later is represented by the consumer's copy
        (above); a product that is NOT consumed becomes its own leaf (the target,
        or a reagent-derived byproduct root).
      * Result: each molecule node has at most one incoming and one outgoing
        reaction edge, so no single molecule "goes to" two reactions.

    Molecule ``role``: "start" (external input, no producer), "target"
    (final product), "byproduct" (reagent-derived root) or "intermediate"
    (produced by one reaction and fed into another). Reaction-node colour is a
    frontend concern (yellow normally, red when changed in the diff view).
    """
    warnings = []
    agents = set()      # canonical SMILES that ever appear as a reagent/agent

    # --- pass 1: clean each reaction; learn what is produced/consumed --------
    per_rxn = []        # index -> (reactants, products) | None
    meta = []           # index -> (rid, number, mapped, source_smiles)
    produced = set()
    consumed = set()
    produced_by = {}    # canonical SMILES -> list of reaction indices that make it
    for index, rxn in enumerate(synthesis):
        rid = str(rxn.get("reaction_id", index))
        number = index + 1
        mapped = (rxn.get("full_reaction_smiles_mapped") or "").strip()
        source_smiles = mapped or rxn.get("full_reaction_smiles") or rxn.get("reaction_smiles") or ""
        meta.append((rid, number, mapped, source_smiles))

        sparts = split_rxn_smiles(source_smiles) if source_smiles else []
        if len(sparts) >= 3:
            for a in sparts[1].split("."):
                ca = safe_canonicalize(a)
                if ca:
                    agents.add(ca)

        cleaned = None
        if source_smiles:
            try:
                cleaned = _clean_reaction_for_nodes(source_smiles)
            except Exception as e:  # pragma: no cover - defensive
                warnings.append(f"Rxn {number}: failed to parse ({e}).")
        per_rxn.append(cleaned)
        if cleaned is not None:
            reactants, products = cleaned
            consumed.update(reactants)
            produced.update(products)
            for p in products:
                produced_by.setdefault(p, []).append(index)

    # --- pass 2: expand the reaction DAG into a TREE from each root product ---
    # Starting at each final product, walk backwards (product <- reaction <-
    # reactants) making a fresh copy of every node. A building block used by
    # several reactions therefore gets its whole upstream sub-tree duplicated, so
    # each path is an independent chain (A→B→C, A→B→D, A→B→E) instead of one node
    # fanning out to several reactions.
    reactants_of = [cl[0] if cl else [] for cl in per_rxn]
    # per-reaction atom-mapped fragments, for orienting each reactant to the
    # product it forms (keyed by canonical SMILES -> fragment mol with maps).
    # Skipped when align=False (connectivity scan doesn't depict anything).
    rxn_mapped = ([_mapped_frags(meta[i][2]) for i in range(len(synthesis))]
                  if align else [None] * len(synthesis))
    rxn_template = {}
    for index, rxn in enumerate(synthesis):
        rid, number, mapped, _ = meta[index]
        rxn_template[index] = {
            "type": "reaction", "index": index, "reaction_id": rid,
            "number": number, "label": f"Rxn {number}",
            "raw_description": rxn.get("raw_description", "") or "",
            "reagents": reagent_summary(rxn),
            "has_mapped": bool(mapped),
            "mapping_confidence": rxn.get("mapping_confidence"),
            "invalid": False,
        }

    nodes = []
    edges = []
    used_rxn = set()          # reaction indices that appear in the tree
    counter = [0]
    NODE_CAP = 2500           # guard against pathologically convergent blow-up
    capped = [False]

    def _nid(prefix):
        counter[0] += 1
        return f"{prefix}::{counter[0]}"

    def expand(smiles, path, is_root, node_mol, wedge_in=None):
        """Emit a molecule node for `smiles` and, if some reaction not already on
        this branch produces it, that reaction + its reactants (recursively).

        `node_mol` is this molecule already laid out (oriented to the product it
        forms). Its reactants are then oriented to IT, via the producing
        reaction's atom mapping, so a shared scaffold keeps one orientation down
        the whole tree."""
        prods = [j for j in produced_by.get(smiles, []) if j not in path]
        if is_root:
            role = "byproduct" if smiles in agents else "target"
        elif prods:
            role = "intermediate"
        else:
            role = "start"
        mid = _nid("m")
        node = {"id": mid, "type": "molecule", "smiles": smiles, "role": role}
        eff_wedges = []
        if node_mol is not None:
            node["molblock"] = Chem.MolToMolBlock(node_mol)
            # Consistent wedge/dash: default per stereocentre, overridden by the
            # bond inherited from the parent so a shared centre is drawn the same
            # way. Rendered on demand (apply_wedges); stored as atom-index pairs.
            eff_wedges = _wedge_plan(node_mol, wedge_in or [])
            if eff_wedges:
                node["wedges"] = eff_wedges
            node["n_atoms"] = node_mol.GetNumAtoms()   # heavy-atom count (fallback sizing)
            # 2D bounding box of the (horizontalized) layout, so the node can be a
            # rectangle matching the molecule's shape at a constant bond scale.
            conf = node_mol.GetConformer()
            xs = [conf.GetAtomPosition(i).x for i in range(node_mol.GetNumAtoms())]
            ys = [conf.GetAtomPosition(i).y for i in range(node_mol.GetNumAtoms())]
            node["ext_w"] = round(max(xs) - min(xs), 3) if xs else 0.0
            node["ext_h"] = round(max(ys) - min(ys), 3) if ys else 0.0
            # Compact (abbreviated) depiction: its own extent (for node sizing).
            # CondenseMolAbbreviations KEEPS the aligned conformer of the surviving
            # atoms, so we do NOT recompute coords — the compact depiction stays in
            # the same orientation as the full structure (so a click, which renders
            # the full molblock, matches the tree). The image is condensed on demand
            # at render time (labels don't survive a MolBlock roundtrip).
            amol, albls = abbreviate_mol(node_mol)
            if amol is not None:
                ac = amol.GetConformer()
                axs = [ac.GetAtomPosition(i).x for i in range(amol.GetNumAtoms())]
                ays = [ac.GetAtomPosition(i).y for i in range(amol.GetNumAtoms())]
                node["ext_w_a"] = round(max(axs) - min(axs), 3) if axs else 0.0
                node["ext_h_a"] = round(max(ays) - min(ays), 3) if ays else 0.0
                node["abbr_labels"] = albls
        nodes.append(node)
        if not prods or len(nodes) > NODE_CAP:
            capped[0] = capped[0] or len(nodes) > NODE_CAP
            return mid
        producer = max(prods)                         # the reaction that made it
        rxid = _nid("rx")
        nodes.append({**rxn_template[producer], "id": rxid})
        used_rxn.add(producer)
        edges.append({"source": rxid, "target": mid})           # reaction -> product
        # This molecule's atoms keyed by the producer's atom-map numbers, so its
        # reactants can be oriented to the matching atoms.
        rm = rxn_mapped[producer]
        map2idx = (_map_to_idx(node_mol, rm["prod"][smiles])
                   if node_mol is not None and rm and smiles in rm["prod"] else None)
        idx2map = {v: k for k, v in map2idx.items()} if map2idx else {}
        # A reactant repeated within the same reaction (e.g. a dimerization
        # A.A>>B, or 4x MeI) is ONE branch with `mult` arrows into the reaction —
        # not `mult` duplicated sub-trees.
        for r, mult in Counter(reactants_of[producer]).items():
            child_mol = None
            child_wedge = None
            if align:
                # A small starting material (leaf, <10 heavy atoms) reads better in
                # its own clean layout than squeezed into the product's orientation,
                # so skip alignment for it.
                r_is_leaf = not any(j not in (path | {producer})
                                    for j in produced_by.get(r, []))
                r_mol = Chem.MolFromSmiles(r)
                r_small = r_mol is not None and r_mol.GetNumHeavyAtoms() < 10
                if r_mol is not None and _is_bridged(r_mol):
                    child_mol = _bridged_coords(r)    # cage: 2D-from-3D, not aligned
                elif r_is_leaf and r_small:
                    child_mol = _plain_coords(r)
                elif map2idx is not None and rm and r in rm["react"]:
                    child_mol = _aligned_child(r, rm["react"][r], map2idx, node_mol)
                    # carry this molecule's wedge choices onto the reactant, matched
                    # by atom-map, so a shared stereocentre is drawn identically.
                    cmap = _map_to_idx(child_mol, rm["react"][r]) if child_mol is not None else None
                    if cmap and eff_wedges and idx2map:
                        child_wedge = [(cmap[idx2map[sc]], cmap[idx2map[nbr]])
                                       for sc, nbr in eff_wedges
                                       if idx2map.get(sc) in cmap and idx2map.get(nbr) in cmap]
                if child_mol is None:                 # no usable mapping -> plain layout
                    child_mol = _plain_coords(r)
            cid = expand(r, path | {producer}, False, child_mol, child_wedge)
            for _ in range(mult):
                edges.append({"source": cid, "target": rxid})   # reactant -> reaction
        return mid

    for root in sorted(p for p in produced if p not in consumed):
        root_mol = None
        if align:
            m = Chem.MolFromSmiles(root)
            root_mol = (_bridged_coords(root) if m is not None and _is_bridged(m)
                        else None) or _plain_coords(root)
        expand(root, set(), True, root_mol)
    if capped[0]:
        warnings.append(f"Tree truncated at {NODE_CAP} nodes (highly convergent route).")

    # reactions never reached above (salt/complex, unparseable, or orphaned) —
    # still show them so nothing silently disappears.
    orphan_reactions = []   # cleanable reactions that couldn't be linked into a tree
    for index in range(len(synthesis)):
        if index in used_rxn:
            continue
        rid, number, mapped, source_smiles = meta[index]
        node = {**rxn_template[index], "id": f"rxn::{index}"}
        if per_rxn[index] is None:
            if source_smiles and check_if_reactant_product_same(source_smiles):
                warnings.append(f"Rxn {number}: reactants == products (salt/complex); shown with no molecule edges.")
            elif not source_smiles:
                warnings.append(f"Rxn {number}: no reaction SMILES available.")
            else:
                node["invalid"] = True
                warnings.append(f"Rxn {number}: SMILES could not be parsed — shown as invalid (red).")
        else:
            orphan_reactions.append(number)   # valid reaction, product/reactant mismatch
        nodes.append(node)

    n_targets = sum(1 for n in nodes
                    if n["type"] == "molecule" and n["role"] == "target")
    # If nothing is a terminal product it's usually a cycle (a step remakes a
    # molecule consumed upstream); surface it so "0 targets" is explainable.
    cycles = _reaction_cycles(per_rxn) if n_targets == 0 else []
    return {
        "nodes": nodes,
        "edges": edges,
        "warnings": warnings,
        "n_steps": len(synthesis),
        "n_targets": n_targets,
        "orphan_reactions": orphan_reactions,
        "connected": n_targets == 1 and not orphan_reactions,
        "cycles": cycles,
    }
