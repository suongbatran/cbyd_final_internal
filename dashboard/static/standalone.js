"use strict";

// STANDALONE build of the read-only viewer (see dashboard/build_static.py).
// Same UI as app.js, but with no server: route data is embedded in the page and
// every structure is drawn in the browser by RDKit.js. Adds "Report a mistake"
// links (e-mail / GitHub issue / pull request) configured by window.DASH_CONFIG.

cytoscape.use(cytoscapeDagre);

const CONFIG = Object.assign(
  { repo: "", branch: "main", email: "", dataPath: "data/routes", rdkitWasmUrl: "", built: "" },
  window.DASH_CONFIG || {});

const ROLE_COLOR = {
  start: "#2e9e5b",
  intermediate: "#2f7fd1",
  target: "#e0533d",
  byproduct: "#9aa0ac",
};
const RXN_YELLOW = "#f4c430";

// Reagent/condition categories. The first four carry a paired *_DEFINITIONS
// SMILES list (DEF_FOR); the rest are name-only.
const REAGENT_FIELDS = [
  ["GENERAL_REAGENTS", "Reagents"],
  ["CATALYSTS", "Catalysts"],
  ["LIGANDS", "Ligands"],
  ["SOLVENTS", "Solvents"],
  ["MODIFIERS", "Modifiers"],
  ["TEMPERATURE", "Temperature"],
  ["YIELD", "Yield"],
];
const DEF_FOR = {
  GENERAL_REAGENTS: "GENERAL_REAGENTS_DEFINITIONS",
  CATALYSTS: "CATALYSTS_DEFINITIONS",
  LIGANDS: "LIGANDS_DEFINITIONS",
  SOLVENTS: "SOLVENTS_DEFINITIONS",
};
const BG_FIELDS = [
  ["name", "name"], ["author", "author"], ["year", "year"],
  ["drug name", "drug name"], ["manufacturer", "manufacturer"],
  ["citation", "citation"], ["doi", "doi"],
];

let RDKit = null;
let CY = null;
let CURRENT = null;
let LAST_MOL = null;    // last molecule node shown (so a toggle can refresh the panel)
let PANEL = null;       // {kind:'reaction'|'molecule'|'info'|'report', ...} — for refresh

// ---- embedded data ---------------------------------------------------------
const INDEX = JSON.parse(document.getElementById("route-index").textContent);
const ROUTE_BLOBS = JSON.parse(document.getElementById("route-data").textContent);
const BY_NAME = new Map(INDEX.map((r) => [r.name, r]));

function b64ToBytes(b64) {
  const bin = atob(b64);
  const out = new Uint8Array(bin.length);
  for (let i = 0; i < bin.length; i++) out[i] = bin.charCodeAt(i);
  return out;
}
async function gunzipJson(b64) {
  const stream = new Blob([b64ToBytes(b64)]).stream()
    .pipeThrough(new DecompressionStream("gzip"));
  return new Response(stream).json();
}

// Unpack one route into the same shape the server's /api/routes/{name} returned.
async function loadRoute(name) {
  const d = await gunzipJson(ROUTE_BLOBS[name]);
  const reactions = d.rx.map((r, i) => ({
    index: i, number: i + 1, reaction_id: r.rid,
    raw_description: r.d, full_reaction_smiles_mapped: r.m, parsed_reaction: r.p || {},
  }));
  const nodes = d.nodes.map((n) => n.t === "m"
    ? { id: n.id, type: "molecule", smiles: n.s, role: n.r, cx: n.cx, cxa: n.cxa, wd: !!n.wd,
        ext_w: n.ew, ext_h: n.eh, ext_w_a: n.ewa, ext_h_a: n.eha, n_atoms: n.na }
    : { id: n.id, type: "reaction", index: n.i, label: `Rxn ${n.i + 1}`,
        raw_description: reactions[n.i].raw_description, invalid: !!n.inv });
  return {
    name,
    background: d.bg || {},
    reactions,
    graph: { nodes, edges: d.edges.map(([s, t]) => ({ source: s, target: t })), warnings: d.warn || [] },
  };
}

// ---- drawing (RDKit.js) ----------------------------------------------------
// Whether to draw heteroatoms in their element colours (O red, N blue, …). Default
// off = monochrome (all black).
let COLOR = false;
// Whether to condense groups/chains (Bn, Boc, PEG6, (CH₂)ₙ…). Default on.
let ABBREV = true;
const BW_PALETTE = (() => { const p = {}; for (let z = -1; z <= 118; z++) p[z] = [0, 0, 0]; return p; })();
const paletteOpts = () => (COLOR ? {} : { atomColourPalette: BW_PALETTE });

const IMG_CACHE = new Map();
const IMG_CACHE_MAX = 6000;
function cached(key, make) {
  let v = IMG_CACHE.get(key);
  if (v === undefined) {
    v = make();
    if (IMG_CACHE.size >= IMG_CACHE_MAX) IMG_CACHE.delete(IMG_CACHE.keys().next().value);
    IMG_CACHE.set(key, v);
  }
  return v;
}
const svgUri = (svg) => "data:image/svg+xml;charset=utf-8," + encodeURIComponent(svg);
function placeholderSvg(text, w, h) {
  return `<svg xmlns="http://www.w3.org/2000/svg" width="${w}" height="${h}" viewBox="0 0 ${w} ${h}">` +
    `<rect x="0.5" y="0.5" width="${w - 1}" height="${h - 1}" fill="#f5f5f5" stroke="#c8c8c8"/>` +
    `<text x="50%" y="50%" fill="#969696" font-family="sans-serif" font-size="${Math.max(10, Math.round(h / 12))}" ` +
    `text-anchor="middle" dominant-baseline="middle">${escapeHtml(text)}</text></svg>`;
}

// Draw a molecule from SMILES / CXSMILES. A CXSMILES from the build carries the
// tree-aligned 2D coords, so RDKit.js keeps that orientation. `wedged` = honour
// the stored wU/wD wedges instead of recomputing them (as the server did).
function molSvg(input, w, h, wedged = false) {
  return cached(`m|${input}|${w}|${h}|${COLOR}|${wedged}`, () => {
    let mol = null;
    try {
      mol = input ? RDKit.get_mol(input) : null;
      if (!mol || !mol.is_valid()) return placeholderSvg("invalid", w, h);
      const details = Object.assign({ width: w, height: h, padding: 0.08 }, paletteOpts(),
        wedged ? { useMolBlockWedging: true, addChiralHs: false } : {});
      return mol.get_svg_with_highlights(JSON.stringify(details));
    } catch (e) {
      return placeholderSvg("invalid", w, h);
    } finally {
      if (mol) mol.delete();
    }
  });
}
const molUrl = (smiles, size = 200) => svgUri(molSvg(smiles, size, size));
// Tree node depiction: aligned coords; the abbreviated form when ABBREV is on.
const nodeUrl = (d, w, h, abbr) =>
  d.cx ? svgUri(molSvg(abbr && d.cxa ? d.cxa : d.cx, w, h, d.wd)) : molUrl(d.smiles, Math.max(w, h));

const useAbbr = (d) => ABBREV && d && d.ext_w_a != null && !!d.cxa;
const NODE_SS = Math.min(3, Math.max(2, Math.ceil(window.devicePixelRatio || 1) + 1));
const NODE_MIN_PX = 260;
const nodeBgImg = (d) => {
  let w = (d.w || 76) * NODE_SS, h = (d.h || 76) * NODE_SS;
  const long = Math.max(w, h);
  if (long < NODE_MIN_PX) { const k = NODE_MIN_PX / long; w *= k; h *= k; }
  return nodeUrl(d, Math.round(w), Math.round(h), useAbbr(d));
};
// Big, aspect-preserving render for the panel / lightbox. Always the FULL structure.
const mbBig = (d, long) => {
  const ew = d.ext_w || 1, eh = d.ext_h || 1;
  if (d.cx && (ew || eh)) {
    return ew >= eh ? nodeUrl(d, long, Math.max(60, Math.round(long * eh / ew)), false)
                    : nodeUrl(d, Math.max(60, Math.round(long * ew / eh)), long, false);
  }
  return nodeUrl(d, long, long, false);
};

// Mapped reaction image (same sizing rule as the server's _reaction_png).
function rxnSvg(smiles, sub = 320) {
  const subW = Math.round(sub * 1.25), subH = sub;
  if (!smiles) return placeholderSvg("no reaction smiles", 700, 240);
  return cached(`r|${smiles}|${sub}|${COLOR}`, () => {
    let rxn = null;
    try {
      const parts = smiles.split(/(?<!-)>(?!-)/);
      const nMain = parts.length >= 2
        ? parts[0].split(".").filter(Boolean).length + parts[parts.length - 1].split(".").filter(Boolean).length
        : 2;
      const width = Math.min(4000, Math.max(2, nMain + 1) * subW);
      rxn = RDKit.get_rxn(smiles, JSON.stringify({ useSmiles: true }));
      if (!rxn) throw new Error("could not parse reaction");
      // CoordGen blanks some reactions (same workaround as the server).
      if (RDKit.prefer_coordgen) RDKit.prefer_coordgen(false);
      try {
        return rxn.get_svg_with_highlights(JSON.stringify(Object.assign({ width, height: subH }, paletteOpts())));
      } finally {
        if (RDKit.prefer_coordgen) RDKit.prefer_coordgen(true);
      }
    } catch (e) {
      return placeholderSvg("invalid reaction", Math.max(subW, 400), subH);
    } finally {
      if (rxn) rxn.delete();
    }
  });
}
const rxnUrl = (smiles, sub = 320) => svgUri(rxnSvg(smiles, sub));
let _fullSizeBlobUrl = null;
function rxnBlobUrl(smiles, sub) {
  if (_fullSizeBlobUrl) URL.revokeObjectURL(_fullSizeBlobUrl);
  _fullSizeBlobUrl = URL.createObjectURL(new Blob([rxnSvg(smiles, sub)], { type: "image/svg+xml" }));
  return _fullSizeBlobUrl;
}

// ---- route list ------------------------------------------------------------
const _FILTER_TEXT_IDS = ["route-bgname", "route-doi", "route-author", "route-year"];
function saveFilters() {
  const f = {};
  for (const id of _FILTER_TEXT_IDS) f[id] = document.getElementById(id).value;
  try { sessionStorage.setItem("dash-view-filters", JSON.stringify(f)); } catch (e) {}
}
function restoreFilters() {
  let f;
  try { f = JSON.parse(sessionStorage.getItem("dash-view-filters") || "{}"); } catch (e) { f = {}; }
  for (const id of _FILTER_TEXT_IDS) if (f[id] != null) document.getElementById(id).value = f[id];
}
const val = (id) => document.getElementById(id).value.toLowerCase().trim();
const low = (v) => (v == null ? "" : String(v).toLowerCase());

// Compound / DOI / author / year — each an ANDed, case-insensitive substring filter.
function loadRouteList() {
  const bg = val("route-bgname"), doi = val("route-doi"), au = val("route-author"), yr = val("route-year");
  const items = INDEX.filter((r) =>
    (!bg || r.q_compound.includes(bg)) && (!doi || r.q_doi.includes(doi)) &&
    (!au || low(r.author).includes(au)) && (!yr || low(r.year).includes(yr)));
  const ul = document.getElementById("route-list");
  const frag = document.createDocumentFragment();
  document.getElementById("route-count").textContent = `${items.length} shown of ${INDEX.length} routes`;
  for (const r of items) {
    const li = document.createElement("li");
    li.dataset.name = r.name;
    const cname = (r.title && r.title !== r.name) ? r.title : null;
    const sub = [cname, r.author, r.year, r.n_steps != null ? `${r.n_steps} rxns` : null]
      .filter(Boolean).join(" · ");
    li.innerHTML =
      `<div class="r-title"><span class="r-num">#${r.number}</span> ${escapeHtml(r.name)}.json</div>` +
      `<div class="r-sub">${escapeHtml(sub)}</div>`;
    li.addEventListener("click", () => selectRoute(r.name));
    frag.appendChild(li);
  }
  ul.innerHTML = "";
  ul.appendChild(frag);
  if (CURRENT) markActive(CURRENT.name);
}
function markActive(name) {
  document.querySelectorAll("#route-list li").forEach((li) =>
    li.classList.toggle("active", li.dataset.name === name));
}

// ---- route load ------------------------------------------------------------
let _routeSeq = 0;
async function selectRoute(name) {
  if (!ROUTE_BLOBS[name]) return;
  const seq = ++_routeSeq;
  markActive(name);
  const route = await loadRoute(name);
  if (seq !== _routeSeq) return;   // a newer click won
  CURRENT = route;
  try { history.replaceState(null, "", "#route=" + encodeURIComponent(name)); } catch (e) {}
  renderHeader();
  showSingle();
  clearPanel();
}

function renderHeader() {
  renderAbbrevUI();
  renderColorUI();
  const bg = CURRENT.background || {};
  document.getElementById("route-title").textContent = `${CURRENT.name}.json`;
  document.getElementById("route-meta").textContent =
    [bg.name, bg.author, bg.year].filter(Boolean).join(" · ");
  document.getElementById("route-tools").classList.remove("hidden");

  const warnEl = document.getElementById("warnings");
  const warns = (CURRENT.graph && CURRENT.graph.warnings) || [];
  if (warns.length) {
    warnEl.innerHTML = `<b>${warns.length} note(s):</b><ul>` +
      warns.map((w) => `<li>${escapeHtml(w)}</li>`).join("") + "</ul>";
    warnEl.classList.remove("hidden");
  } else warnEl.classList.add("hidden");
}

// ---- graph rendering -------------------------------------------------------
function graphStyle() {
  return [
    {
      selector: 'node[type="molecule"]',
      style: {
        "shape": "round-rectangle",
        "width": (ele) => ele.data("w") || 76,
        "height": (ele) => ele.data("h") || 76,
        "background-color": "#fff",
        "background-image": (ele) => nodeBgImg(ele.data()),
        "background-fit": "contain", "background-clip": "node",
        "border-width": 3, "border-color": (ele) => ROLE_COLOR[ele.data("role")] || "#999",
      },
    },
    {
      selector: 'node[type="reaction"]',
      style: {
        "shape": "round-rectangle",
        "label": (ele) => (ele.data("invalid") ? "⚠ " : "") + ele.data("label"),
        "width": 62, "height": 30,
        "background-color": (ele) => ele.data("invalid") ? "#ffd9d2" : RXN_YELLOW,
        "color": "#222",
        "font-size": 11, "font-weight": 700,
        "text-valign": "center", "text-halign": "center",
        "border-width": (ele) => ele.data("invalid") ? 3 : 2,
        "border-style": (ele) => ele.data("invalid") ? "dashed" : "solid",
        "border-color": (ele) => ele.data("invalid") ? "#e0533d" : "#c79a14",
      },
    },
    { selector: "node:selected", style: { "border-color": "#111", "border-width": 4 } },
    {
      selector: "edge",
      style: {
        "width": 2, "line-color": "#b9bdc7", "target-arrow-color": "#b9bdc7",
        "target-arrow-shape": "triangle", "curve-style": "bezier", "arrow-scale": 0.9,
      },
    },
    {
      selector: 'edge[kind="rxn2rxn"]',
      style: { "line-color": "#8a8f9c", "target-arrow-color": "#8a8f9c", "width": 3 },
    },
  ];
}

function molNodeDims(n) {
  const SCALE = 16, MIN = 46, DEFAULT = 76;
  const abbr = useAbbr(n);
  const ew = abbr ? n.ext_w_a : n.ext_w, eh = abbr ? n.ext_h_a : n.ext_h;
  if ((ew || eh) && ew != null) {
    return { w: Math.max(MIN, Math.round(ew * SCALE)),
             h: Math.max(MIN, Math.round(eh * SCALE)) };
  }
  const s = n.n_atoms ? Math.max(MIN, Math.round(19 * Math.sqrt(n.n_atoms))) : DEFAULT;
  return { w: s, h: s };
}

function toElements(graph) {
  const els = [];
  const rxnIds = new Set(graph.nodes.filter(n => n.type === "reaction").map(n => n.id));
  for (const n of graph.nodes) {
    if (n.type === "molecule") {
      const dim = molNodeDims(n);
      els.push({ data: { id: n.id, type: "molecule", smiles: n.smiles, role: n.role,
                         cx: n.cx, cxa: n.cxa, wd: n.wd, w: dim.w, h: dim.h,
                         ext_w: n.ext_w, ext_h: n.ext_h, ext_w_a: n.ext_w_a, ext_h_a: n.ext_h_a } });
    } else {
      els.push({ data: {
        id: n.id, type: "reaction", index: n.index, label: n.label,
        raw: n.raw_description, invalid: !!n.invalid,
      }});
    }
  }
  graph.edges.forEach((e, i) => {
    const kind = (rxnIds.has(e.source) && rxnIds.has(e.target)) ? "rxn2rxn" : "mol";
    els.push({ data: { id: `e${i}__${e.source}__${e.target}`, source: e.source, target: e.target, kind } });
  });
  return els;
}

function bindWheelScroll(container) {
  if (container.dataset.wheelBound) return;
  container.dataset.wheelBound = "1";
  container.addEventListener("wheel", (e) => {
    const cy = container.__cy;
    if (!cy) return;
    e.preventDefault();
    if (e.ctrlKey || e.metaKey) {
      const level = Math.min(Math.max(cy.zoom() * Math.pow(1.01, -e.deltaY),
                                      cy.minZoom()), cy.maxZoom());
      const r = container.getBoundingClientRect();
      cy.zoom({ level, renderedPosition: { x: e.clientX - r.left, y: e.clientY - r.top } });
    } else if (e.shiftKey) {
      cy.panBy({ x: -(e.deltaY || e.deltaX), y: 0 });
    } else {
      cy.panBy({ x: -e.deltaX, y: -e.deltaY });
    }
  }, { passive: false });
}

function spaceComponents(cy) {
  const comps = cy.elements().components();
  if (comps.length <= 1) return;
  const items = comps.map((c) => ({ c, bb: c.boundingBox() })).sort((a, b) => a.bb.x1 - b.bb.x1);
  const maxW = Math.max(...items.map((it) => it.bb.w));
  const gap = Math.min(Math.max(maxW * 0.7, 260), 900);
  let right = null;
  for (const { c, bb } of items) {
    if (right === null) { right = bb.x2; continue; }
    const dx = (right + gap) - bb.x1;
    if (dx) c.nodes().positions((n) => ({ x: n.position("x") + dx, y: n.position("y") }));
    right = bb.x2 + dx;
  }
}

function makeGraph(containerId, graph, fit = true) {
  const container = document.getElementById(containerId);
  const cy = cytoscape({
    container,
    elements: toElements(graph),
    minZoom: 0.06, maxZoom: 3,
    userZoomingEnabled: false,
    style: graphStyle(),
    layout: { name: "dagre", rankDir: "BT", nodeSep: 24, rankSep: 28, edgeSep: 10, fit: false, padding: 20 },
  });
  spaceComponents(cy);
  if (fit) cy.fit(undefined, 20);
  container.__cy = cy;
  bindWheelScroll(container);
  cy.on("tap", 'node[type="reaction"]', (evt) => showReaction(evt.target.data("index")));
  cy.on("tap", 'node[type="molecule"]', (evt) => showMolecule(evt.target.data()));
  cy.on("mouseover", "node", (evt) => {
    const d = evt.target.data();
    if (d.type === "molecule") { hideTooltip(); return; }
    showTooltip(d.raw || d.label, evt.renderedPosition, containerId);
  });
  cy.on("mousemove", 'node[type="reaction"]', (evt) => moveTooltip(evt.renderedPosition, containerId));
  cy.on("mouseout", "node", hideTooltip);
  return cy;
}

function renderAbbrevUI() {
  const btn = document.getElementById("btn-abbr");
  if (btn) btn.textContent = "Abbrev: " + (ABBREV ? "on" : "off");
}
function renderColorUI() {
  const btn = document.getElementById("btn-color");
  if (btn) btn.textContent = "Colored atoms: " + (COLOR ? "on" : "off");
}
// Re-render whichever panel is open, so its image reflects a toggle.
function refreshPanel() {
  if (!PANEL) return;
  if (PANEL.kind === "reaction") showReaction(PANEL.index);
  else if (PANEL.kind === "molecule" && LAST_MOL) showMolecule(LAST_MOL);
  else if (PANEL.kind === "info") showRouteInfo();
}
function rerenderGraph() {
  if (CURRENT) showSingle();
}

function showSingle() {
  if (CY) CY.destroy();
  CY = makeGraph("cy", CURRENT.graph, true);
}

// ---- panel -----------------------------------------------------------------
function clearPanel() {
  PANEL = null;
  document.getElementById("panel-empty").classList.remove("hidden");
  document.getElementById("panel-content").classList.add("hidden");
}
function openPanel(html) {
  const el = document.getElementById("panel-content");
  el.innerHTML = html;
  document.getElementById("panel-empty").classList.add("hidden");
  el.classList.remove("hidden");
  bindReportBox(el);
  return el;
}

function showMolecule(d) {
  LAST_MOL = d;
  PANEL = { kind: "molecule" };
  const role = d.role;
  openPanel(
    `<div class="panel-step" style="color:${ROLE_COLOR[role] || "#666"}">${role}</div>` +
    `<img class="panel-img zoomable" src="${mbBig(d, 520)}" alt="molecule" title="click to enlarge" />` +
    `<div class="section-title">SMILES</div>` +
    `<div class="smiles-box" style="min-height:auto">${escapeHtml(d.smiles || "∅")}</div>` +
    reportBox({ kind: "molecule", role, smiles: d.smiles }, "Report a mistake in this molecule"));
  document.querySelector("#panel-content .panel-img")
    .addEventListener("click", () => openLightbox(mbBig(d, 2000)));
}

function showReaction(index) {
  const rxn = CURRENT.reactions.find((r) => r.index === index);
  if (!rxn) return;
  PANEL = { kind: "reaction", index };
  const mapped = rxn.full_reaction_smiles_mapped || "";
  const parsed = rxn.parsed_reaction || {};

  // Reagents/conditions: name card + a "↳ SMILES" card with structure thumbnails.
  // Empty categories are skipped.
  let reagentSection = "";
  for (const [key, label] of REAGENT_FIELDS) {
    const names = parsed[key] || [];
    const defKey = DEF_FOR[key];
    const list = defKey ? (parsed[defKey] || []) : [];
    if (!names.length && !list.length) continue;
    let card = `<div class="reagent-card"><div class="fld"><span class="fld-label">${label}</span>` +
      `<div class="smiles-box" style="min-height:auto">${escapeHtml(names.join("; ")) || "∅"}</div></div>`;
    if (defKey) {
      const thumbs = list.map((s) =>
        `<img class="reagent-thumb" src="${molUrl(s, 80)}" data-mol="${escapeHtml(s)}" title="${escapeHtml(s)}" alt="${escapeHtml(s)}"/>`).join("");
      card += `<div class="fld fld-smiles"><span class="fld-label">↳ SMILES</span>` +
        (thumbs ? `<div class="thumb-row">${thumbs}</div>` : "") +
        `<div class="smiles-box" style="min-height:auto">${escapeHtml(list.join("; ")) || "∅"}</div></div>`;
    }
    reagentSection += card + `</div>`;
  }
  if (!reagentSection) reagentSection = `<div class="muted">No reagent information.</div>`;

  const roBox = (label, v) =>
    `<div class="fld"><span class="fld-label">${label}</span>` +
    `<div class="smiles-box">${escapeHtml(v || "∅")}</div></div>`;

  openPanel(
    `<div class="panel-step">Rxn ${rxn.number} <span class="muted">· id ${escapeHtml(String(rxn.reaction_id))} · index ${rxn.index}</span></div>` +
    `<img class="panel-img zoomable" id="rxn-img" src="${rxnUrl(mapped, 360)}" alt="reaction" title="click to enlarge" />` +
    `<div class="img-actions">` +
      `<a id="rxn-full" href="#" target="_blank" rel="noopener">Open full size ↗</a>` +
    `</div>` +
    `<div class="section-title">Description</div>` +
    `<div class="smiles-box" style="min-height:auto">${escapeHtml(rxn.raw_description || "∅")}</div>` +
    `<div class="section-title">Reagents &amp; conditions</div>` + reagentSection +
    `<div class="section-title">Reaction SMILES</div>` +
    roBox("full_reaction_smiles_mapped", mapped) +
    reportBox({ kind: "reaction", rxn }, "Report a mistake in this reaction"));

  document.getElementById("rxn-img").addEventListener("click", () => openLightbox(rxnUrl(mapped, 820)));
  document.getElementById("rxn-full").addEventListener("click", (e) => {
    e.currentTarget.href = rxnBlobUrl(mapped, 760);
  });
}

// ---- panel: route info (read-only background) ------------------------------
function showRouteInfo() {
  if (!CURRENT) return;
  PANEL = { kind: "info" };
  const bg = CURRENT.background || {};
  let rows = "";
  for (const [key, label] of BG_FIELDS) {
    const v = bg[key] == null ? "" : String(bg[key]);
    rows += `<div class="fld"><span class="fld-label">${label}</span>` +
      `<div class="smiles-box" style="min-height:auto">${escapeHtml(v || "∅")}</div></div>`;
  }
  openPanel(
    `<div class="panel-step">Route info — ${escapeHtml(CURRENT.name)}.json</div>` + rows +
    reportBox({ kind: "route" }, "Report a mistake in this route"));
}

function showReport() {
  if (!CURRENT) return;
  PANEL = { kind: "report" };
  openPanel(
    `<div class="panel-step">Report a mistake — ${escapeHtml(CURRENT.name)}.json</div>` +
    `<div class="panel-desc">Spotted a wrong structure, reagent, condition, or reference? ` +
    `Tip: click the reaction or molecule first, then use the report box in its panel ` +
    `so the details are filled in for you.</div>` +
    reportBox({ kind: "route" }, null));
}

// ---- reporting (e-mail / GitHub issue / pull request) ----------------------
// Every link is pre-filled with the route, file, and (if any) reaction/molecule
// details plus the reporter's note. GitHub's file editor turns an edit into a PR
// (forking first for people without write access).
const REPORTS = new Map();
let _reportSeq = 0;

const routeFilePath = () => `${CONFIG.dataPath}/${CURRENT.name}.json`;
const ghFileUrl = (mode) =>
  `https://github.com/${CONFIG.repo}/${mode}/${encodeURIComponent(CONFIG.branch)}/` +
  routeFilePath().split("/").map(encodeURIComponent).join("/");

function reportText(scope, note) {
  const bg = CURRENT.background || {};
  const num = (BY_NAME.get(CURRENT.name) || {}).number;
  let subject = `[CbyD] Mistake in ${CURRENT.name}`;
  const lines = [
    "**What is wrong / suggested fix:**",
    note.trim() || "_(please describe the problem and, if you can, the correct value)_",
    "",
    "---",
    `- Route: \`${CURRENT.name}.json\`` + (num ? ` (#${num})` : ""),
    `- File: \`${routeFilePath()}\``,
  ];
  if (bg.name) lines.push(`- Compound: ${bg.name}`);
  const ref = [bg.author, bg.year, bg.citation, bg.doi ? `doi:${bg.doi}` : null].filter(Boolean).join(", ");
  if (ref) lines.push(`- Reference: ${ref}`);
  if (scope.kind === "reaction") {
    const r = scope.rxn;
    subject += ` · Rxn ${r.number}`;
    lines.push(`- Reaction: Rxn ${r.number} (reaction_id ${r.reaction_id}, index ${r.index})`);
    if (r.raw_description) lines.push(`- Description: ${r.raw_description}`);
    if (r.full_reaction_smiles_mapped) lines.push(`- Mapped SMILES: \`${r.full_reaction_smiles_mapped}\``);
  } else if (scope.kind === "molecule") {
    subject += ` · ${scope.role} molecule`;
    lines.push(`- Molecule (${scope.role}): \`${scope.smiles}\``);
  }
  if (CONFIG.built) lines.push(`- Viewer build: ${CONFIG.built}`);
  return { subject, body: lines.join("\n") };
}

function reportLinks(scope, note) {
  const { subject, body } = reportText(scope, note);
  const enc = encodeURIComponent;
  const links = {};
  if (CONFIG.email) {
    links.email = `mailto:${CONFIG.email}?subject=${enc(subject)}&body=${enc(body)}`;
  }
  if (CONFIG.repo) {
    // Browsers / GitHub cap URL length; trim the body rather than fail.
    let b = body;
    const base = `https://github.com/${CONFIG.repo}/issues/new?title=${enc(subject)}&body=`;
    while (base.length + enc(b).length > 7500 && b.length > 200) b = b.slice(0, Math.floor(b.length * 0.8)) + "…";
    links.issue = base + enc(b);
    links.pr = ghFileUrl("edit");
    links.view = ghFileUrl("blob");
  }
  return links;
}

function reportBox(scope, summary) {
  const id = `report-${++_reportSeq}`;
  REPORTS.set(id, scope);
  const configured = CONFIG.email || CONFIG.repo;
  const inner = configured
    ? `<textarea placeholder="What's wrong, and what should it be? (optional — you can also write it in the e-mail / issue)"></textarea>` +
      `<div class="report-actions">` +
        (CONFIG.email ? `<a class="btn btn-primary" data-link="email" href="#">✉ Email</a>` : "") +
        (CONFIG.repo ? `<a class="btn btn-primary" data-link="issue" href="#" target="_blank" rel="noopener">Open GitHub issue ↗</a>` +
                       `<a class="btn" data-link="pr" href="#" target="_blank" rel="noopener" ` +
                       `title="Edit ${escapeHtml(CURRENT.name)}.json on GitHub and open a pull request">Suggest a fix (PR) ↗</a>` +
                       `<a class="btn" data-link="view" href="#" target="_blank" rel="noopener">View JSON ↗</a>` : "") +
      `</div>` +
      (CONFIG.repo ? `<div class="report-ctx" style="margin-top:6px">PR: GitHub opens <code>${escapeHtml(routeFilePath())}</code> ` +
                     `in its editor — make the change, then “Propose changes”.</div>` : "")
    : `<div class="muted">Reporting links are not configured for this copy of the viewer.</div>`;
  const box = `<div class="report-box" data-report="${id}">${inner}</div>`;
  return summary
    ? `<details class="report-details" style="margin-top:16px"><summary class="section-title" style="cursor:pointer">${escapeHtml(summary)}</summary>${box}</details>`
    : box;
}

function bindReportBox(root) {
  root.querySelectorAll(".report-box[data-report]").forEach((box) => {
    const scope = REPORTS.get(box.dataset.report);
    const ta = box.querySelector("textarea");
    const update = () => {
      const links = reportLinks(scope, ta ? ta.value : "");
      box.querySelectorAll("a[data-link]").forEach((a) => { a.href = links[a.dataset.link] || "#"; });
    };
    if (ta) ta.addEventListener("input", update);
    update();
  });
}

// ---- lightbox --------------------------------------------------------------
function openLightbox(src) {
  document.getElementById("lightbox-img").src = src;
  document.getElementById("lightbox").classList.remove("hidden");
}
function closeLightbox() {
  document.getElementById("lightbox").classList.add("hidden");
  document.getElementById("lightbox-img").src = "";
}

// ---- tooltip ---------------------------------------------------------------
const tip = document.getElementById("tooltip");
function showTooltip(content, pos, containerId) {
  if (!content) return;
  tip.textContent = content;
  tip.classList.remove("img-tip");
  tip.classList.remove("hidden");
  moveTooltip(pos, containerId);
}
function moveTooltip(pos, containerId) {
  const rect = document.getElementById(containerId || "cy").getBoundingClientRect();
  tip.style.left = rect.left + pos.x + 14 + "px";
  tip.style.top = rect.top + pos.y + 14 + "px";
}
function hideTooltip() { tip.classList.add("hidden"); }

// Hover a reagent thumbnail (right panel) -> larger structure preview near cursor.
function _placeMolHover(x, y) {
  const w = 312;
  let left = x + 18, top = y + 18;
  if (left + w > window.innerWidth) left = x - w - 18;
  if (top + w > window.innerHeight) top = Math.max(8, window.innerHeight - w - 8);
  tip.style.left = left + "px";
  tip.style.top = top + "px";
}
const _panelEl = document.getElementById("panel");
if (_panelEl) {
  _panelEl.addEventListener("mouseover", (e) => {
    const t = e.target.closest && e.target.closest(".reagent-thumb");
    if (!t || !t.dataset.mol) return;
    tip.innerHTML = `<img src="${molUrl(t.dataset.mol, 320)}" alt="reagent" />`;
    tip.classList.add("img-tip");
    tip.classList.remove("hidden");
    _placeMolHover(e.clientX, e.clientY);
  });
  _panelEl.addEventListener("mousemove", (e) => {
    if (tip.classList.contains("hidden")) return;
    if (e.target.closest && e.target.closest(".reagent-thumb")) _placeMolHover(e.clientX, e.clientY);
  });
  _panelEl.addEventListener("mouseout", (e) => {
    if (e.target.closest && e.target.closest(".reagent-thumb")) hideTooltip();
  });
}

// ---- utils -----------------------------------------------------------------
function escapeHtml(s) {
  return String(s == null ? "" : s)
    .replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;").replace(/"/g, "&quot;");
}
function debounce(fn, ms) { let t; return (...a) => { clearTimeout(t); t = setTimeout(() => fn(...a), ms); }; }

// ---- resizable columns -----------------------------------------------------
function resizeGraphs() { if (CY) CY.resize(); }
function initResizers() {
  const app = document.getElementById("app");
  const getVar = (n, d) => parseFloat(getComputedStyle(app).getPropertyValue(n)) || d;
  try {
    const s = JSON.parse(localStorage.getItem("dash-cols") || "{}");
    if (s.left) app.style.setProperty("--left", s.left + "px");
    if (s.right) app.style.setProperty("--right", s.right + "px");
  } catch (e) {}

  let drag = null;
  document.querySelectorAll(".gutter").forEach((g) => {
    g.addEventListener("mousedown", (e) => {
      drag = { side: g.dataset.side, x: e.clientX, el: g,
               left: getVar("--left", 280), right: getVar("--right", 460) };
      g.classList.add("dragging");
      document.body.style.cursor = "col-resize";
      document.body.style.userSelect = "none";
      e.preventDefault();
    });
  });
  window.addEventListener("mousemove", (e) => {
    if (!drag) return;
    const dx = e.clientX - drag.x;
    const maxW = window.innerWidth - 360;
    if (drag.side === "left") {
      app.style.setProperty("--left", Math.min(Math.max(drag.left + dx, 180), maxW) + "px");
    } else {
      app.style.setProperty("--right", Math.min(Math.max(drag.right - dx, 240), maxW) + "px");
    }
    resizeGraphs();
  });
  window.addEventListener("mouseup", () => {
    if (!drag) return;
    drag.el.classList.remove("dragging");
    document.body.style.cursor = "";
    document.body.style.userSelect = "";
    try {
      localStorage.setItem("dash-cols",
        JSON.stringify({ left: getVar("--left", 280), right: getVar("--right", 460) }));
    } catch (e) {}
    drag = null;
    resizeGraphs();
  });
}

// ---- init ------------------------------------------------------------------
async function initRDKit() {
  const wasmEl = document.getElementById("rdkit-wasm");
  const b64 = wasmEl ? wasmEl.textContent.trim() : "";
  const opts = b64 ? { wasmBinary: b64ToBytes(b64) } : { locateFile: () => CONFIG.rdkitWasmUrl };
  if (wasmEl) wasmEl.remove();   // free the base64 text
  RDKit = await window.initRDKitModule(opts);
  if (RDKit.prefer_coordgen) RDKit.prefer_coordgen(true);   // cleaner 2D for SMILES-only drawings
}

function renderAbout() {
  const bits = [];
  if (CONFIG.built) bits.push(`Built ${escapeHtml(CONFIG.built)}`);
  if (CONFIG.repo) bits.push(`<a href="https://github.com/${escapeHtml(CONFIG.repo)}" target="_blank" rel="noopener">GitHub</a>`);
  document.getElementById("about").innerHTML = bits.join(" · ");
}

const _reloadList = debounce(() => { saveFilters(); loadRouteList(); }, 200);
for (const id of _FILTER_TEXT_IDS) document.getElementById(id).addEventListener("input", _reloadList);
document.getElementById("btn-info").addEventListener("click", () => CURRENT && showRouteInfo());
document.getElementById("btn-report").addEventListener("click", showReport);
document.getElementById("btn-abbr").addEventListener("click", () => {
  ABBREV = !ABBREV;
  renderAbbrevUI();
  rerenderGraph();
});
document.getElementById("btn-color").addEventListener("click", () => {
  COLOR = !COLOR;
  renderColorUI();
  rerenderGraph();
  refreshPanel();
});
document.getElementById("lightbox").addEventListener("click", closeLightbox);
document.addEventListener("keydown", (e) => { if (e.key === "Escape") closeLightbox(); });

(async function main() {
  initResizers();
  restoreFilters();
  renderAbout();
  loadRouteList();
  try {
    await initRDKit();
  } catch (e) {
    const el = document.getElementById("loading");
    el.classList.add("error");
    document.getElementById("loading-msg").textContent =
      "Could not start the chemistry engine (RDKit.js): " + e.message +
      (CONFIG.rdkitWasmUrl ? " — this copy loads it from the internet; check your connection." : "");
    return;
  }
  document.getElementById("loading").remove();
  // Deep link: #route=<name> opens that route (handy in bug reports).
  const m = location.hash.match(/route=([^&]+)/);
  if (m) {
    const name = decodeURIComponent(m[1]);
    if (ROUTE_BLOBS[name]) selectRoute(name);
  }
})();
