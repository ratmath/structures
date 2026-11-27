// Groups given by generators and relations, and Todd–Coxeter coset
// enumeration.
//
// A presentation ⟨a, b | a⁴, b², (ab)²⟩ names generators and words that are
// to equal the identity (relators). The group is the free group on the
// generators modulo everything those words force. Coset enumeration builds
// the action of the generators on the right cosets of a subgroup H, given by
// words; when it finishes, the number of cosets is the index [G : H], and with
// H = 1 that is the order of G and the action is the group itself.
//
// Words are arrays of nonzero integers: k stands for generator k (from 1) and
// −k for its inverse.

import { Group, sup } from "./group.js";
import { fromRelations } from "./abelian.js";

const SUPER = { "⁰": "0", "¹": "1", "²": "2", "³": "3", "⁴": "4", "⁵": "5", "⁶": "6", "⁷": "7", "⁸": "8", "⁹": "9", "⁻": "-" };

/** Free reduction: cancel x x⁻¹ pairs. */
export function reduceWord(w) {
  const out = [];
  for (const x of w) {
    if (out.length && out[out.length - 1] === -x) out.pop();
    else out.push(x);
  }
  return out;
}
export const invertWord = (w) => w.map((x) => -x).reverse();

function cyclicallyReduce(w) {
  let r = reduceWord(w);
  while (r.length > 1 && r[0] === -r[r.length - 1]) r = r.slice(1, -1);
  return r;
}

/** A word as text: a²b⁻¹ (or a^2*b^-1 with { plain: true }), "e" when empty. */
export function formatWord(w, names, { plain = false } = {}) {
  if (!w.length) return "e";
  const parts = [];
  let i = 0;
  while (i < w.length) {
    let j = i;
    while (j < w.length && w[j] === w[i]) j += 1;
    const k = (j - i) * Math.sign(w[i]);
    const name = names[Math.abs(w[i]) - 1];
    parts.push(k === 1 ? name : plain ? `${name}^${k}` : `${name}${sup(k)}`);
    i = j;
  }
  return parts.join(plain ? "*" : "");
}

/** Reads a word over the given generator names: "a^2 b^-1", "a²b⁻¹", "(ab)^3", "[a, b]" (= a⁻¹b⁻¹ab), "A" for a⁻¹. */
export function parseWord(text, names) {
  const src = String(text).replace(/[⁰¹²³⁴⁵⁶⁷⁸⁹⁻]+/g, (m) => `^${m.split("").map((c) => SUPER[c]).join("")}`).replace(/−/g, "-");
  const sorted = names.map((n, i) => [n, i]).sort((a, b) => b[0].length - a[0].length);
  let pos = 0;
  const peek = () => src[pos];
  const skip = () => { while (pos < src.length && /[\s*·]/.test(src[pos])) pos += 1; };
  const power = () => {
    skip();
    if (peek() !== "^") return 1;
    pos += 1;
    skip();
    const m = /^\(?\s*(-?\d+)\s*\)?/.exec(src.slice(pos));
    if (!m) throw new SyntaxError(`a whole-number power is needed at "${src.slice(pos)}"`);
    pos += m[0].length;
    return Number(m[1]);
  };
  const pw = (w, k) => {
    const base = k < 0 ? invertWord(w) : w;
    const out = [];
    for (let i = 0; i < Math.abs(k); i += 1) out.push(...base);
    return out;
  };
  const atom = () => {
    skip();
    const c = peek();
    if (c === "(") {
      pos += 1;
      const w = seq();
      skip();
      if (peek() !== ")") throw new SyntaxError("a ( is not closed");
      pos += 1;
      return w;
    }
    if (c === "[") {
      pos += 1;
      const a = seq();
      skip();
      if (peek() !== ",") throw new SyntaxError("a commutator [x, y] needs a comma");
      pos += 1;
      const b = seq();
      skip();
      if (peek() !== "]") throw new SyntaxError("a [ is not closed");
      pos += 1;
      return [...invertWord(a), ...invertWord(b), ...a, ...b];
    }
    if (c === "1" || c === "e" && !names.includes("e")) {
      if (c === "1" || !/[A-Za-z0-9_]/.test(src[pos + 1] ?? "")) { pos += 1; return []; }
    }
    for (const [n, i] of sorted) if (src.startsWith(n, pos)) { pos += n.length; return [i + 1]; }
    // a capital letter for the inverse of a one-letter lowercase generator
    if (c && c !== c.toLowerCase()) {
      const i = names.indexOf(c.toLowerCase());
      if (i >= 0) { pos += 1; return [-(i + 1)]; }
    }
    throw new SyntaxError(`unknown generator at "${src.slice(pos)}" (generators: ${names.join(", ")})`);
  };
  const seq = () => {
    const out = [];
    for (;;) {
      skip();
      const c = peek();
      if (c === undefined || c === ")" || c === "]" || c === "," || c === "=") return out;
      const w = atom();
      out.push(...pw(w, power()));
    }
  };
  const w = seq();
  skip();
  if (pos < src.length) throw new SyntaxError(`could not read "${src.slice(pos)}"`);
  return reduceWord(w);
}

const splitTop = (s, sep) => {
  const out = [];
  let depth = 0;
  let cur = "";
  for (const ch of s) {
    if (ch === "(" || ch === "[") depth += 1;
    if (ch === ")" || ch === "]") depth -= 1;
    if (ch === sep && depth === 0) { out.push(cur); cur = ""; continue; }
    cur += ch;
  }
  out.push(cur);
  return out.map((x) => x.trim()).filter(Boolean);
};

export class Presentation {
  /** Generators by name and relators as words or text ("a^4", "(ab)^2", "ab = ba"). */
  constructor(generators, relators = []) {
    if (new Set(generators).size !== generators.length) throw new RangeError("two generators have the same name");
    this.generators = Object.freeze([...generators]);
    this.relators = Object.freeze(relators.map((r) => Object.freeze(this.#relator(r))).filter((w) => w.length));
    Object.freeze(this);
  }

  #relator(r) {
    if (Array.isArray(r)) return cyclicallyReduce(r);
    const sides = String(r).split("=");
    if (sides.length > 2) throw new SyntaxError(`${r}: one = at most`);
    const lhs = parseWord(sides[0], this.generators);
    const rhs = sides.length === 2 ? parseWord(sides[1], this.generators) : [];
    return cyclicallyReduce([...lhs, ...invertWord(rhs)]);
  }

  /** A word from text, over these generators. */
  word(text) { return Array.isArray(text) ? reduceWord(text) : parseWord(text, this.generators); }
  formatWord(w, options) { return formatWord(w, this.generators, options); }

  /**
   * Todd–Coxeter coset enumeration (HLT strategy, with coincidences) of the
   * right cosets of the subgroup generated by the words H. Throws when more
   * than maxCosets cosets are defined, which is what happens for an infinite
   * index.
   */
  enumerate(H = [], { maxCosets = 500000 } = {}) {
    return enumerate(this.generators.length, this.relators, H.map((h) => this.word(h)), maxCosets, this.generators);
  }

  /** [G : H]. */
  index(H = [], options) { return this.enumerate(H, options).index; }
  /** |G|, when coset enumeration over the trivial subgroup finishes. */
  order(options) { return this.enumerate([], options).index; }

  /**
   * The group as a Group (it must be finite): elements are the permutations
   * of the cosets of 1, printed as shortest words in the generators.
   */
  toGroup({ maxCosets, name = null } = {}) {
    const T = this.enumerate([], { maxCosets });
    if (T.index > 20000) throw new RangeError(`order ${T.index} is too large to list`);
    const perms = this.generators.map((_, g) => T.table.map((row) => row[g]));
    const reps = T.representatives;
    const compose = (a, b) => a.map((x) => b[x]); // right action: first a, then b
    const id = reps.map((_, i) => i);
    const G = Group.generate({
      generators: perms,
      identity: id,
      mul: compose,
      key: (p) => (Array.isArray(p) ? p.join(",") : String(p)),
      format: (p) => this.formatWord(reps[p[0]]),
      name: name ?? this.toString(),
      labels: Object.fromEntries(reps.map((w, i) => [this.formatWord(w, { plain: true }), T.elementPermutation(i)])),
    });
    G.presentation = this;
    G.fromWord = (w) => T.elementPermutation(T.cosetOfWord(this.word(w)));
    return G;
  }

  /**
   * The abelianisation G/G′ as an AbelianGroup: the relators, counted by
   * exponent sums, are relations of a Z-module. Works for infinite groups too.
   */
  abelianInvariants() {
    const n = this.generators.length;
    const rows = this.relators.map((w) => {
      const v = Array(n).fill(0n);
      for (const x of w) v[Math.abs(x) - 1] += x > 0 ? 1n : -1n;
      return v;
    });
    return fromRelations(rows, n);
  }

  toString() {
    return `⟨${this.generators.join(", ")} | ${this.relators.map((w) => this.formatWord(w)).join(", ")}⟩`;
  }
}

/** A presentation from text "a, b | a^4, b^2, (ab)^2" (angle brackets optional), or from names and relators. */
export function presentation(gens, relators) {
  if (relators === undefined && typeof gens === "string") {
    const t = gens.trim().replace(/^[⟨<]\s*/, "").replace(/\s*[⟩>]$/, "");
    const [g, r = ""] = t.split("|");
    const names = g.split(",").map((s) => s.trim()).filter(Boolean);
    return new Presentation(names, splitTop(r, ","));
  }
  return new Presentation(typeof gens === "string" ? gens.split(",").map((s) => s.trim()) : gens, relators ?? []);
}

// ---- Todd–Coxeter -----------------------------------------------------------------------------------

function enumerate(k, relators, subgroup, maxCosets, names) {
  const cols = 2 * k; // column 2i: generator i; column 2i+1: its inverse
  const colOf = (x) => (x > 0 ? 2 * (x - 1) : 2 * (-x - 1) + 1);
  const rels = relators.map((w) => w.map(colOf));
  const subs = subgroup.map((w) => w.map(colOf));
  let cap = 1024;
  let T = new Int32Array(cap * cols).fill(-1);
  let parent = new Int32Array(cap);
  let n = 1; // cosets defined so far (live or dead)
  parent[0] = 0;
  let total = 1;
  const grow = () => {
    const c2 = cap * 2;
    const T2 = new Int32Array(c2 * cols).fill(-1);
    T2.set(T);
    const p2 = new Int32Array(c2);
    p2.set(parent);
    T = T2;
    parent = p2;
    cap = c2;
  };
  const define = (c, x) => {
    if (total >= maxCosets) throw new RangeError(`coset enumeration stopped after ${maxCosets} cosets: the index may be infinite, or larger than maxCosets`);
    if (n >= cap) grow();
    const d = n;
    n += 1;
    total += 1;
    parent[d] = d;
    T[d * cols + x] = -1;
    for (let j = 0; j < cols; j += 1) T[d * cols + j] = -1;
    T[c * cols + x] = d;
    T[d * cols + (x ^ 1)] = c;
  };
  const rep = (c) => {
    let r = c;
    while (parent[r] !== r) r = parent[r];
    while (parent[c] !== r) { const nx = parent[c]; parent[c] = r; c = nx; }
    return r;
  };
  const queue = [];
  const merge = (a, b) => {
    a = rep(a); b = rep(b);
    if (a === b) return;
    if (a > b) [a, b] = [b, a];
    parent[b] = a;
    queue.push(b);
  };
  const coincidence = (a, b) => {
    merge(a, b);
    while (queue.length) {
      const g = queue.shift();
      for (let x = 0; x < cols; x += 1) {
        const d = T[g * cols + x];
        if (d < 0) continue;
        T[d * cols + (x ^ 1)] = -1;
        const mu = rep(g);
        const nu = rep(d);
        if (T[mu * cols + x] >= 0) merge(nu, T[mu * cols + x]);
        else if (T[nu * cols + (x ^ 1)] >= 0) merge(mu, T[nu * cols + (x ^ 1)]);
        else { T[mu * cols + x] = nu; T[nu * cols + (x ^ 1)] = mu; }
      }
      total -= 1;
    }
  };
  const live = (c) => parent[c] === c;
  const scanAndFill = (c, w) => {
    const r = w.length;
    let f = c;
    let i = 0;
    let b = c;
    let j = r - 1;
    for (;;) {
      while (i <= j && T[f * cols + w[i]] >= 0) { f = T[f * cols + w[i]]; i += 1; }
      if (i > j) { if (f !== b) coincidence(f, b); return; }
      while (j >= i && T[b * cols + (w[j] ^ 1)] >= 0) { b = T[b * cols + (w[j] ^ 1)]; j -= 1; }
      if (j < i) { coincidence(f, b); return; }
      if (i === j) { T[f * cols + w[i]] = b; T[b * cols + (w[i] ^ 1)] = f; return; }
      define(f, w[i]);
    }
  };
  for (const h of subs) if (h.length) scanAndFill(0, h);
  for (let a = 0; a < n; a += 1) {
    for (const w of rels) {
      if (!live(a)) break;
      scanAndFill(a, w);
    }
    if (!live(a)) continue;
    for (let x = 0; x < cols; x += 1) if (T[a * cols + x] < 0) define(a, x);
  }
  // renumber the live cosets in breadth-first order from coset 0, which also gives shortest words
  const order = new Int32Array(n).fill(-1);
  const list = [0];
  order[0] = 0;
  const repWord = [[]];
  for (let i = 0; i < list.length; i += 1) {
    const c = list[i];
    for (let x = 0; x < cols; x += 1) {
      const d = T[c * cols + x];
      if (d < 0 || order[d] >= 0) continue;
      order[d] = list.length;
      list.push(d);
      repWord.push([...repWord[i], x % 2 ? -(x >> 1) - 1 : (x >> 1) + 1]);
    }
  }
  const index = list.length;
  const table = list.map((c) => Array.from({ length: k }, (_, g) => order[T[c * cols + 2 * g]]));
  const inverseTable = list.map((c) => Array.from({ length: k }, (_, g) => order[T[c * cols + 2 * g + 1]]));
  for (const row of table) if (row.some((x) => x < 0)) throw new Error("internal: coset table not complete");
  const cosetOfWord = (w) => {
    let c = 0;
    for (const x of w) c = x > 0 ? table[c][x - 1] : inverseTable[c][-x - 1];
    return c;
  };
  return Object.freeze({
    index,
    table,
    inverseTable,
    representatives: repWord,
    defined: n,
    cosetOfWord,
    /** With H = 1: the permutation of cosets by which the element of coset i acts (on the right) — i.e. coset j ↦ j·w_i. */
    elementPermutation: (i) => table.map((_, j) => walk(j, repWord[i])),
    toString: () => `${index} cosets; generators act by\n${names.map((nm, g) => `${nm}: ${table.map((r) => r[g]).join(" ")}`).join("\n")}`,
  });
  function walk(c, w) {
    for (const x of w) c = x > 0 ? table[c][x - 1] : inverseTable[c][-x - 1];
    return c;
  }
}

/**
 * A presentation of a finite Group: its generators, and one relator for each
 * edge of its Cayley graph outside a spanning tree (Schreier's method). Long
 * but always correct; coset enumeration on it gives back the group.
 */
export function presentationOf(G, { names = null } = {}) {
  const k = G.gens.length;
  const nm = names ?? (k <= 26 ? Array.from({ length: k }, (_, i) => String.fromCharCode(97 + i)) : Array.from({ length: k }, (_, i) => `g${i + 1}`));
  const { order, from, by } = G.spanningTree();
  const word = new Array(G.order);
  word[0] = [];
  for (const x of order) if (x !== 0) word[x] = [...word[from[x]], by[x] + 1];
  const rels = [];
  const seen = new Set();
  for (let x = 0; x < G.order; x += 1) for (let s = 0; s < k; s += 1) {
    const y = G.mulIndex(x, G.gens[s]);
    if (from[y] === x && by[y] === s) continue; // a tree edge
    const w = cyclicallyReduce([...word[x], s + 1, ...invertWord(word[y])]);
    if (!w.length) continue;
    const key = w.join(",");
    if (seen.has(key)) continue;
    seen.add(key);
    rels.push(w);
  }
  return new Presentation(nm, rels);
}
