"use strict";

// Read-only VIEWER build of the dashboard. Same tree/graph rendering as the edit
// build, but every editing/review/staging control is removed — you can only browse
// routes and inspect reactions & molecules.

cytoscape.use(cytoscapeDagre);

const ROLE_COLOR = {
  start: "#2e9e5b",
  intermediate: "#2f7fd1",
  target: "#e0533d",
  byproduct: "#9aa0ac",
};
const RXN_YELLOW = "#f4c430";

// Reagent/condition categories (same layout as the edit build). The first four
// carry a paired *_DEFINITIONS SMILES list (DEF_FOR); the rest are name-only.
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

let CY = null;
let CURRENT = null;
let LAST_MOL = null;    // last molecule node shown (so a toggle can refresh the panel)
let PANEL = null;       // {kind:'reaction'|'molecule'|'info', ...} — current panel for refresh

// ---- API -------------------------------------------------------------------
async function api(path, opts) {
  const res = await fetch(path, opts);
  if (!res.ok) {
    let detail = res.statusText;
    try { detail = (await res.json()).detail || detail; } catch (e) {}
    throw new Error(detail);
  }
  return res.json();
}

// Whether to draw heteroatoms in their element colours (O red, N blue, …). Default
// off = monochrome (all black). Appended to every image URL so colour/mono are
// cached separately server-side and the browser refetches when it changes.
let COLOR = false;
const colorParam = () => (COLOR ? "&color=1" : "");
const molUrl = (smiles, size = 200) =>
  `/api/mol.png?smiles=${encodeURIComponent(smiles)}&size=${size}` + colorParam();
const mbUrl = (key, w = 200, h = 0, abbr = false) =>
  `/api/mol_mb.png?key=${encodeURIComponent(key)}&w=${w}&h=${h || w}` +
  (abbr ? "&abbr=1" : "") + colorParam();
// Whether to condense groups/chains (Bn, Boc, PEG6, (CH₂)ₙ…). Default on.
let ABBREV = true;
const useAbbr = (d) => ABBREV && d && d.ext_w_a != null;
const NODE_SS = Math.min(3, Math.max(2, Math.ceil(window.devicePixelRatio || 1) + 1));
const NODE_MIN_PX = 260;
const nodeBgImg = (d) => {
  let w = (d.w || 76) * NODE_SS, h = (d.h || 76) * NODE_SS;
  const long = Math.max(w, h);
  if (long < NODE_MIN_PX) { const k = NODE_MIN_PX / long; w *= k; h *= k; }
  return d && d.mb ? mbUrl(d.mb, Math.round(w), Math.round(h), useAbbr(d))
                   : molUrl(d.smiles, Math.round(Math.max(w, h)));
};
// Big, aspect-preserving render for the panel / lightbox. Always the FULL structure.
const mbBig = (d, long) => {
  const ew = d.ext_w || 1, eh = d.ext_h || 1;
  if (d.mb && (ew || eh)) {
    return ew >= eh ? mbUrl(d.mb, long, Math.max(60, Math.round(long * eh / ew)))
                    : mbUrl(d.mb, Math.max(60, Math.round(long * ew / eh)), long);
  }
  return d.mb ? mbUrl(d.mb, long, long) : molUrl(d.smiles, long);
};
const rxnUrl = (smiles, sub = 320) =>
  `/api/reaction.png?smiles=${encodeURIComponent(smiles)}&sub=${sub}` + colorParam();

// ---- route list ------------------------------------------------------------
function currentQuery() {
  return {
    bg: document.getElementById("route-bgname").value,
    doi: document.getElementById("route-doi").value,
    author: document.getElementById("route-author").value,
    year: document.getElementById("route-year").value,
  };
}
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
let _listSeq = 0;
async function loadRouteList() {
  const { bg, doi, author, year } = currentQuery();
  const seq = ++_listSeq;
  const data = await api(
    `/api/routes?bg=${encodeURIComponent(bg || "")}&doi=${encodeURIComponent(doi || "")}` +
    `&author=${encodeURIComponent(author || "")}&year=${encodeURIComponent(year || "")}&limit=0`);
  if (seq !== _listSeq) return;   // a newer filter won
  const ul = document.getElementById("route-list");
  ul.innerHTML = "";
  document.getElementById("route-count").textContent = `${data.count} shown of ${data.n_all} routes`;
  for (const r of data.routes) {
    const li = document.createElement("li");
    li.dataset.name = r.name;
    const cname = (r.title && r.title !== r.name) ? r.title : null;
    const sub = [cname, r.author, r.year, r.n_steps != null ? `${r.n_steps} rxns` : null]
      .filter(Boolean).join(" · ");
    li.innerHTML =
      `<div class="r-title"><span class="r-num">#${r.number}</span> ${escapeHtml(r.name)}.json</div>` +
      `<div class="r-sub">${escapeHtml(sub)}</div>`;
    li.addEventListener("click", () => selectRoute(r.name));
    ul.appendChild(li);
  }
}
function markActive(name) {
  document.querySelectorAll("#route-list li").forEach((li) =>
    li.classList.toggle("active", li.dataset.name === name));
}

// ---- route load ------------------------------------------------------------
async function selectRoute(name) {
  markActive(name);
  CURRENT = await api(`/api/routes/${encodeURIComponent(name)}`);
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
  const abbr = ABBREV && n.ext_w_a != null;
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
      els.push({ data: { id: n.id, type: "molecule", smiles: n.smiles, role: n.role, mb: n.mb,
                         w: dim.w, h: dim.h,
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
    `<div class="smiles-box" style="min-height:auto">${escapeHtml(d.smiles || "∅")}</div>`);
  document.querySelector("#panel-content .panel-img")
    .addEventListener("click", () => openLightbox(mbBig(d, 2000)));
}

function showReaction(index) {
  const rxn = CURRENT.reactions.find((r) => r.index === index);
  if (!rxn) return;
  PANEL = { kind: "reaction", index };
  const mapped = rxn.full_reaction_smiles_mapped || "";
  const parsed = rxn.parsed_reaction || {};

  // Reagents/conditions: same card layout as the edit build (name card + a "↳ SMILES"
  // card with structure thumbnails), but read-only. Empty categories are skipped.
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

  const roBox = (label, val) =>
    `<div class="fld"><span class="fld-label">${label}</span>` +
    `<div class="smiles-box">${escapeHtml(val || "∅")}</div></div>`;

  const bigUrl = rxnUrl(mapped, 820);
  openPanel(
    `<div class="panel-step">Rxn ${rxn.number} <span class="muted">· id ${escapeHtml(String(rxn.reaction_id))} · index ${rxn.index}</span></div>` +
    `<img class="panel-img zoomable" id="rxn-img" src="${rxnUrl(mapped, 360)}" alt="reaction" title="click to enlarge" />` +
    `<div class="img-actions">` +
      `<a href="${rxnUrl(mapped, 760)}" target="_blank" rel="noopener">Open full size ↗</a>` +
    `</div>` +
    `<div class="section-title">Description</div>` +
    `<div class="smiles-box" style="min-height:auto">${escapeHtml(rxn.raw_description || "∅")}</div>` +
    `<div class="section-title">Reagents &amp; conditions</div>` + reagentSection +
    `<div class="section-title">Reaction SMILES</div>` +
    roBox("full_reaction_smiles_mapped", mapped));

  document.getElementById("rxn-img").addEventListener("click", () => openLightbox(bigUrl));
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
    `<div class="panel-step">Route info — ${escapeHtml(CURRENT.name)}.json</div>` + rows);
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
    localStorage.setItem("dash-cols",
      JSON.stringify({ left: getVar("--left", 280), right: getVar("--right", 460) }));
    drag = null;
    resizeGraphs();
  });
}

// ---- init ------------------------------------------------------------------
const _reloadList = debounce(() => { saveFilters(); loadRouteList(); }, 200);
document.getElementById("route-bgname").addEventListener("input", _reloadList);
document.getElementById("route-doi").addEventListener("input", _reloadList);
document.getElementById("route-author").addEventListener("input", _reloadList);
document.getElementById("route-year").addEventListener("input", _reloadList);
document.getElementById("btn-info").addEventListener("click", () => CURRENT && showRouteInfo());
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

initResizers();
restoreFilters();
loadRouteList();
