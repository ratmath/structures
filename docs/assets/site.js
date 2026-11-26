// Shared page furniture for the structures docs: the contents column, the
// theme switch, result rendering and expression demos (any element with
// data-expr becomes one).

import * as A from "./structures.js?v=f408e7fbe3";
import { makeLang } from "./lang.js?v=f408e7fbe3";
import { mountRunboxes } from "./runbox.js?v=f408e7fbe3";
import { makeGroups } from "./groups.js?v=f408e7fbe3";

const { run } = makeLang(A);
const groups = makeGroups(A);
export { A, run, groups };

export const PAGES = [
  { href: "index.html", title: "Overview" },
  { href: "reference.html", title: "Reference" },
  { href: "explore.html", title: "Explore" },
  { group: "Groups" },
  { href: "groups.html", title: "Groups and how to build them", blurb: "Cyclic, dihedral, symmetric, quaternion and matrix groups; products; Cayley tables" },
  { href: "classes.html", title: "Conjugacy classes and commutators", blurb: "The centre, the class equation, solvable and nilpotent groups" },
  { href: "subgroups.html", title: "Subgroups and quotients", blurb: "Lagrange's theorem, normal subgroups, cosets, quotient groups" },
  { href: "sylow.html", title: "Sylow subgroups", blurb: "Subgroups of prime-power order and how many there are" },
  { href: "homomorphisms.html", title: "Homomorphisms and isomorphism", blurb: "Kernels, the isomorphism theorems, automorphisms, the groups of order up to 20" },
  { href: "actions.html", title: "Actions and counting", blurb: "Orbits, stabilisers, and counting necklaces and cubes up to symmetry" },
  { group: "Polynomials" },
  { href: "polynomials.html", title: "Polynomials in several variables", blurb: "Rings, exact coefficients, reading and printing" },
  { href: "orders.html", title: "Monomial orders and division", blurb: "lex, grlex, grevlex and dividing by several polynomials" },
  { href: "gcd.html", title: "Greatest common divisors", blurb: "Contents, primitive parts, square-free parts" },
  { href: "factoring.html", title: "Factoring", blurb: "Irreducible factors over Q by Kronecker's trick" },
  { group: "Elimination" },
  { href: "resultants.html", title: "Resultants and discriminants", blurb: "Eliminating a variable; where curves meet" },
  { href: "rational-functions.html", title: "Rational functions", blurb: "The field Q(x, y) and matrices with parameters" },
  { group: "Symmetry" },
  { href: "symmetric.html", title: "Symmetric polynomials", blurb: "Elementary symmetric polynomials and Newton's identities" },
];

export const h = (tag, attrs = {}, ...kids) => {
  const el = document.createElement(tag);
  for (const [k, v] of Object.entries(attrs)) {
    if (v === null || v === undefined || v === false) continue;
    if (k === "class") el.className = v;
    else if (k.startsWith("on")) el.addEventListener(k.slice(2), v);
    else el.setAttribute(k, v === true ? "" : v);
  }
  for (const kid of kids.flat(Infinity)) if (kid !== null && kid !== undefined && kid !== false) el.append(kid instanceof Node ? kid : document.createTextNode(String(kid)));
  return el;
};

// ---- chrome -------------------------------------------------------------------

function theme() {
  let saved = null;
  try {
    saved = localStorage.getItem("ratmath-theme");
  } catch {
    /* storage may be blocked */
  }
  const apply = (t) => {
    if (t === "auto") document.documentElement.removeAttribute("data-theme");
    else document.documentElement.setAttribute("data-theme", t);
    try {
      localStorage.setItem("ratmath-theme", t);
    } catch {
      /* ignore */
    }
    for (const b of document.querySelectorAll(".theme button")) b.classList.toggle("on", b.dataset.t === t);
  };
  if (saved) apply(saved);
  return h("div", { class: "theme", role: "group", "aria-label": "Theme" },
    ["auto", "light", "dark"].map((t) => h("button", { "data-t": t, class: (saved ?? "auto") === t ? "on" : "", onclick: () => apply(t) }, t)));
}

function chrome() {
  const here = location.pathname.split("/").pop() || "index.html";
  const aside = document.querySelector("aside");
  if (!aside) return;
  const nav = h("nav", { class: "toc", "aria-label": "Contents" });
  const sel = h("select", { id: "toc-select", "aria-label": "Page", onchange: (e) => (location.href = e.target.value) });
  let group = null;
  for (const p of PAGES) {
    if (p.group) {
      nav.append(h("div", { class: "group" }, p.group));
      group = h("optgroup", { label: p.group });
      sel.append(group);
      continue;
    }
    nav.append(h("a", { href: p.href, class: p.href === here ? "on" : null, "aria-current": p.href === here ? "page" : null }, p.title));
    (group ?? sel).append(h("option", { value: p.href, selected: p.href === here }, p.title));
  }
  aside.append(
    h("header", { class: "brand" }, h("a", { href: "index.html" }, h("div", { class: "mark" }, "RatMath ", h("span", {}, "structures"))), h("small", {}, "Groups and polynomials, exactly")),
    sel, nav, theme(),
  );
}

/** Topic cards for the overview page. */
export function topicCards(el) {
  el.append(...PAGES.filter((p) => p.blurb).map((p) => h("a", { href: p.href }, h("strong", {}, p.title), h("span", {}, p.blurb))));
}

// ---- results --------------------------------------------------------------------

/** One result descriptor (from lang.run) as DOM. */
export function render(r) {
  if (!r) return h("span", { class: "muted" }, "—");
  switch (r.type) {
    case "text": return h("div", { class: "big" }, r.text);
    case "record": return h("div", {}, h("div", { class: "big" }, r.title), h("dl", { class: "kv" }, r.rows.flatMap(([k, x]) => [h("dt", {}, k), h("dd", {}, x)])));
    case "matrix": {
      const tbl = h("table", { class: "data" });
      if (r.head) tbl.append(h("tr", {}, r.head.map((x) => h("th", {}, x))));
      for (const row of r.rows) tbl.append(h("tr", {}, row.map((x) => h("td", {}, x))));
      return h("div", { class: "tablewrap" }, tbl);
    }
    default: return h("pre", {}, JSON.stringify(r));
  }
}

/** Run `work` after the browser has painted a busy state; shows the time it took. */
export function later(out, work) {
  out.classList.add("busy");
  setTimeout(() => {
    const t0 = performance.now();
    let node;
    try {
      node = work();
      out.classList.remove("err");
    } catch (e) {
      node = h("div", { class: "err" }, e.message);
      out.classList.add("err");
    }
    const ms = performance.now() - t0;
    out.replaceChildren(node, h("div", { class: "timing" }, ms < 1 ? "< 1 ms" : `${ms.toFixed(0)} ms`));
    out.classList.remove("busy");
  }, 20);
}

/** An expression demo: input, examples, result. Lines (or ;-separated parts outside parentheses) run in turn. */
export function exprDemo(el, { expr, examples = [], title = null, rows = 1 }) {
  const input = rows > 1
    ? h("textarea", { class: "text", rows, spellcheck: "false", "aria-label": "Expression" }, expr)
    : h("input", { class: "text grow", value: expr, spellcheck: "false", "aria-label": "Expression" });
  const out = h("div", { class: "out" });
  const go = () => later(out, () => {
    const env = {};
    const lines = input.value.split(/\n/).map((s) => s.trim()).filter(Boolean);
    let last = null;
    for (const line of lines) last = run(line, env).result;
    return render(last);
  });
  input.addEventListener("keydown", (e) => {
    if (e.key === "Enter" && (rows === 1 || e.metaKey || e.ctrlKey)) {
      e.preventDefault();
      go();
    }
  });
  el.classList.add("demo");
  el.replaceChildren(...[
    title ? h("div", { class: "demo-title" }, title) : null,
    h("div", { class: "demo-row" }, input, h("button", { onclick: go }, "Run")),
    examples.length ? h("div", { class: "demo-row" }, examples.map((x) => h("button", { class: "chipbtn", onclick: () => { input.value = x.replace(/\\n/g, "\n"); go(); } }, x))) : null,
    out,
  ].filter(Boolean));
  go();
}

/** A group explorer: a group's name in, its invariants (and, if small, its Cayley table) out. */
export function groupDemo(el, { group, examples = [], title = null, table = false }) {
  const input = h("input", { class: "text grow", value: group, spellcheck: "false", "aria-label": "Group" });
  const out = h("div", { class: "out" });
  const go = () => later(out, () => {
    const G = groups.parse(input.value);
    const { rows } = groups.report(G);
    const kids = [h("div", { class: "big" }, String(G)), h("dl", { class: "kv" }, rows.flatMap(([k, x]) => [h("dt", {}, k), h("dd", {}, x)]))];
    if (table && G.order <= 24) kids.push(cayley(G));
    return h("div", {}, kids);
  });
  input.addEventListener("keydown", (e) => {
    if (e.key === "Enter") { e.preventDefault(); go(); }
  });
  el.classList.add("demo");
  el.replaceChildren(...[
    title ? h("div", { class: "demo-title" }, title) : null,
    h("div", { class: "demo-row" }, input, h("button", { onclick: go }, "Show")),
    examples.length ? h("div", { class: "demo-row" }, examples.map((x) => h("button", { class: "chipbtn", onclick: () => { input.value = x; go(); } }, x))) : null,
    out,
  ].filter(Boolean));
  go();
}

/** A Cayley table with one colour per element, so patterns show. */
function cayley(G) {
  const T = G.cayleyTable();
  const n = T.labels.length;
  const colour = (j) => `hsl(${Math.round((360 * j) / n)} 70% 50% / .22)`;
  const tbl = h("table", { class: "data cayley" });
  tbl.append(h("tr", {}, h("th", {}, "·"), T.labels.map((l) => h("th", {}, l))));
  T.table.forEach((row, i) => tbl.append(h("tr", {}, h("th", {}, T.labels[i]), row.map((j) => h("td", { style: `background:${colour(j)}` }, T.labels[j])))));
  return h("div", { class: "tablewrap" }, tbl);
}

// ---- boot -------------------------------------------------------------------------

chrome();
for (const el of document.querySelectorAll("[data-expr]")) {
  exprDemo(el, {
    expr: el.dataset.expr.replace(/\\n/g, "\n"),
    examples: el.dataset.examples ? el.dataset.examples.split("|") : [],
    title: el.dataset.title ?? null,
    rows: el.dataset.rows ? Number(el.dataset.rows) : 1,
  });
}
for (const el of document.querySelectorAll("[data-group]")) {
  groupDemo(el, {
    group: el.dataset.group,
    examples: el.dataset.groups ? el.dataset.groups.split("|") : [],
    title: el.dataset.title ?? null,
    table: "table" in el.dataset,
  });
}
mountRunboxes(document, { scope: A });
