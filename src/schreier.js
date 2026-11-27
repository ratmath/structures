// Permutation groups too large to list: the Schreier–Sims algorithm.
//
// A group of permutations of n points is held by a base (points b₁, b₂, …
// that only the identity fixes all of) and, for each k, the orbit of bₖ under
// the subgroup fixing b₁, …, bₖ₋₁ together with one permutation carrying bₖ
// to each point of that orbit. The order of the group is the product of the
// orbit lengths, and every element is one product of such carriers, which is
// how membership is tested (by "sifting") without listing anything.
//
// Permutations are arrays of images counted from 0, as in @ratmath/algebraic;
// products are a∘b, first b then a. Points given to and returned by the
// methods are counted from 1, as in cycle notation.

import { cycleString, fromCycles, PermGroup } from "@ratmath/algebraic";
import { permutationGroup } from "./named.js";

const mulP = (a, b) => { const n = a.length; const r = new Array(n); for (let i = 0; i < n; i += 1) r[i] = a[b[i]]; return r; };
const invP = (a) => { const n = a.length; const r = new Array(n); for (let i = 0; i < n; i += 1) r[a[i]] = i; return r; };
const idP = (n) => Array.from({ length: n }, (_, i) => i);
const isId = (a) => { for (let i = 0; i < a.length; i += 1) if (a[i] !== i) return false; return true; };
const samePerm = (a, b) => { for (let i = 0; i < a.length; i += 1) if (a[i] !== b[i]) return false; return true; };

// A word is null (the identity), a signed generator number ±(i+1), or a pair
// [a, b] for the product a∘b. Pairs share their parts, so a long word costs
// little to hold; it is written out only on request.
const mulW = (a, b) => (a === null ? b : b === null ? a : [a, b]);
const invW = (a) => (a === null ? null : typeof a === "number" ? -a : { inv: a });

function wordLength(w, memo = new Map()) {
  if (w === null) return 0n;
  if (typeof w === "number") return 1n;
  if (memo.has(w)) return memo.get(w);
  const out = Array.isArray(w) ? wordLength(w[0], memo) + wordLength(w[1], memo) : wordLength(w.inv, memo);
  memo.set(w, out);
  return out;
}

function flattenWord(w, limit) {
  const out = [];
  const stack = [[w, false]];
  while (stack.length) {
    const [x, inverted] = stack.pop();
    if (x === null) continue;
    if (typeof x === "number") {
      out.push(inverted ? -x : x);
      if (out.length > limit) return null;
    } else if (Array.isArray(x)) {
      // (ab)⁻¹ = b⁻¹a⁻¹; the stack gives back what is pushed last first
      if (inverted) stack.push([x[0], true], [x[1], true]);
      else stack.push([x[1], false], [x[0], false]);
    } else stack.push([x.inv, !inverted]);
  }
  // cancel neighbours g g⁻¹
  const tidy = [];
  for (const g of out) { if (tidy.length && tidy[tidy.length - 1] === -g) tidy.pop(); else tidy.push(g); }
  return tidy;
}

/** A small seeded generator (mulberry32), so that random elements can be repeated. */
function seeded(seed) {
  let a = (Number(seed) >>> 0) || 0x9e3779b9;
  return () => {
    a = (a + 0x6d2b79f5) | 0;
    let t = Math.imul(a ^ (a >>> 15), 1 | a);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

function readPerm(g, n) {
  if (typeof g === "string") return fromCycles(g, n);
  if (!Array.isArray(g) || g.length !== n) throw new RangeError(`a permutation of ${n} points is an array of ${n} images or cycle notation`);
  const seen = new Uint8Array(n);
  for (const x of g) {
    if (!Number.isInteger(x) || x < 0 || x >= n || seen[x]) throw new RangeError(`[${g}] is not a permutation of ${n} points`);
    seen[x] = 1;
  }
  return [...g];
}

const degreeOf = (generators) => Math.max(1, ...generators.map((g) => (Array.isArray(g) ? g.length : Math.max(0, ...[...String(g).matchAll(/\d+/g)].map((m) => Number(m[0]))))));

/**
 * A permutation group held by a base and strong generating set.
 */
export class PermChain {
  #levels = [];
  #n;
  #gens;
  #names;
  #words;
  #order = null;

  /**
   * @param generators permutations (arrays counted from 0, or cycle notation counted from 1)
   * @param options { degree, base (points, counted from 1, to put first), names (of the generators, for words), name }
   */
  constructor(generators, { degree = undefined, base = [], names = null, name = null, words = true } = {}) {
    if (generators instanceof PermGroup) { degree = generators.degree; generators = generators.generators; }
    const n = degree ?? degreeOf(generators);
    if (!Number.isInteger(n) || n < 1) throw new RangeError("the degree is a positive whole number");
    this.#n = n;
    this.#gens = generators.map((g) => readPerm(g, n));
    this.#names = names ? [...names] : this.#gens.map((_, i) => `g${i + 1}`);
    if (this.#names.length !== this.#gens.length) throw new RangeError("one name for each generator");
    this.#words = Boolean(words);
    this.name = name;
    for (const b of base) {
      if (!Number.isInteger(b) || b < 1 || b > n) throw new RangeError(`base point ${b} is not one of 1, …, ${n}`);
      if (!this.#levels.some((L) => L.point === b - 1)) this.#newLevel(b - 1);
    }
    this.#gens.forEach((p, i) => { if (!isId(p)) this.#insert(0, { p, w: this.#words ? i + 1 : null }); });
    // drop trailing levels that turned out trivial and were not asked for
    while (this.#levels.length > base.length && this.#levels[this.#levels.length - 1].points.length === 1) this.#levels.pop();
    Object.freeze(this);
  }

  #newLevel(point) {
    const e = { p: idP(this.#n), w: null };
    this.#levels.push({ point, gens: [], rep: new Map([[point, e]]), inv: new Map([[point, e]]), points: [point] });
  }

  #times(a, b) { return { p: mulP(a.p, b.p), w: this.#words ? mulW(a.w, b.w) : null }; }
  #inverse(a) { return { p: invP(a.p), w: this.#words ? invW(a.w) : null }; }

  /** Strips h by the carriers from `level` on: what is left, and the level where it stopped. */
  #sift(h, level) {
    for (let i = level; i < this.#levels.length; i += 1) {
      if (isId(h.p)) return { rest: h, level: i, done: true };
      const L = this.#levels[i];
      const u = L.inv.get(h.p[L.point]);
      if (!u) return { rest: h, level: i, done: false };
      h = this.#times(u, h);
    }
    return { rest: h, level: this.#levels.length, done: isId(h.p) };
  }

  #insert(level, g) {
    if (level === this.#levels.length) this.#newLevel(g.p.findIndex((x, i) => x !== i));
    const L = this.#levels[level];
    L.gens.push(g);
    const queue = L.points.map((pt) => [pt, g]);
    for (let head = 0; head < queue.length; head += 1) {
      const [pt, s] = queue[head];
      const img = s.p[pt];
      const prod = this.#times(s, L.rep.get(pt));
      if (!L.rep.has(img)) {
        L.rep.set(img, prod);
        L.inv.set(img, this.#inverse(prod));
        L.points.push(img);
        for (const t of L.gens) queue.push([img, t]);
      } else {
        // a Schreier generator: it fixes this level's point, so it must lie in the levels below
        const sch = this.#times(L.inv.get(img), prod);
        if (isId(sch.p)) continue;
        const r = this.#sift(sch, level + 1);
        if (!r.done) this.#insert(level + 1, r.rest);
      }
    }
  }

  /** The number of points permuted. */
  get degree() { return this.#n; }
  /** The generators the group was given, as arrays counted from 0. */
  get generators() { return this.#gens.map((g) => [...g]); }
  get generatorNames() { return [...this.#names]; }
  /** The base, counted from 1: only the identity fixes every one of these points. */
  get base() { return this.#levels.map((L) => L.point + 1); }
  /** The lengths of the basic orbits; their product is the order. */
  get orbitLengths() { return this.#levels.map((L) => L.points.length); }
  /** The basic orbits themselves, counted from 1. */
  get basicOrbits() { return this.#levels.map((L) => L.points.map((x) => x + 1).sort((a, b) => a - b)); }
  /** A strong generating set: for each k its members fixing the first k − 1 base points generate that stabiliser. */
  get strongGenerators() {
    const seen = new Set();
    const out = [];
    for (const L of this.#levels) for (const g of L.gens) { const k = g.p.join(","); if (!seen.has(k)) { seen.add(k); out.push([...g.p]); } }
    return out;
  }

  /** The order, as a BigInt. */
  get order() {
    return this.#levels.reduce((o, L) => o * BigInt(L.points.length), 1n);
  }

  get isTrivial() { return this.#levels.every((L) => L.points.length === 1); }

  #read(g) {
    if (g instanceof Object && !Array.isArray(g) && Array.isArray(g.permutation)) g = g.permutation;
    return readPerm(g, this.#n);
  }

  /** Whether the permutation belongs to the group. */
  contains(g) {
    return this.#sift({ p: this.#read(g), w: null }, 0).done;
  }

  /**
   * Sifts a permutation through the chain: { member, level, images, residue }.
   * `images` are the points (from 1) that the base points are sent to, level
   * by level, which name a member uniquely. For a non-member, `level` is
   * where it fell out and `residue` what was left there, in cycle notation.
   */
  sift(g) {
    let h = this.#read(g);
    const images = [];
    for (let i = 0; i < this.#levels.length; i += 1) {
      const L = this.#levels[i];
      const img = h[L.point];
      const u = L.inv.get(img);
      if (!u) return { member: false, level: i + 1, images, residue: cycleString(h), point: img + 1 };
      images.push(img + 1);
      h = mulP(u.p, h);
    }
    if (!isId(h)) return { member: false, level: this.#levels.length + 1, images, residue: cycleString(h), point: null };
    return { member: true, level: this.#levels.length, images, residue: "()", point: null };
  }

  /**
   * A member written in the generators the group was given:
   * { length, word, text, check }. `word` lists signed generator numbers
   * (−k for the inverse of generator k) to be multiplied as written, so the
   * rightmost acts first; `check` multiplies them back. Words from sifting
   * are far from the shortest, and `word` is left null beyond `limit` terms.
   */
  word(g, { limit = 200000 } = {}) {
    if (!this.#words) throw new Error("this group was built with words: false");
    const p = this.#read(g);
    let h = { p, w: null };
    let w = null;
    for (const L of this.#levels) {
      const img = h.p[L.point];
      const u = L.inv.get(img);
      if (!u) throw new RangeError(`${cycleString(p)} is not in the group`);
      h = { p: mulP(u.p, h.p), w: null };
      // g = rep₁ rep₂ … repₖ
      w = mulW(w, L.rep.get(img).w);
    }
    if (!isId(h.p)) throw new RangeError(`${cycleString(p)} is not in the group`);
    const length = wordLength(w);
    const word = length > BigInt(limit) * 4n ? null : flattenWord(w, limit * 4);
    if (!word || word.length > limit) return { length, word: null, text: null, check: null };
    let back = idP(this.#n);
    const inverses = this.#gens.map(invP);
    for (const k of word) back = mulP(back, k > 0 ? this.#gens[k - 1] : inverses[-k - 1]);
    const text = word.length ? collapse(word).map(([k, e]) => this.#names[Math.abs(k) - 1] + (k < 0 ? (e === 1 ? "⁻¹" : `⁻${supers(e)}`) : e === 1 ? "" : supers(e))).join(" ") : "()";
    return { length: BigInt(word.length), word, text, check: samePerm(back, p) };
  }

  /** The product of generator names or signed numbers, as a permutation: the inverse of `word`. */
  evaluate(word) {
    let out = idP(this.#n);
    const inverses = this.#gens.map(invP);
    for (const k of word) {
      let i = k;
      if (typeof k === "string") {
        const neg = /(⁻¹|'|\^-1)$/.test(k);
        const base = k.replace(/(⁻¹|'|\^-1)$/, "");
        const at = this.#names.indexOf(base);
        if (at < 0) throw new RangeError(`no generator named ${base}`);
        i = neg ? -(at + 1) : at + 1;
      }
      if (!Number.isInteger(i) || i === 0 || Math.abs(i) > this.#gens.length) throw new RangeError(`no generator ${k}`);
      out = mulP(out, i > 0 ? this.#gens[i - 1] : inverses[-i - 1]);
    }
    return out;
  }

  /** A member chosen uniformly at random; `seed` or `random` (a function giving numbers in [0, 1)) makes it repeatable. */
  random({ seed = undefined, random = undefined } = {}) {
    const next = random ?? (seed === undefined ? Math.random : seeded(seed));
    let out = idP(this.#n);
    for (const L of this.#levels) out = mulP(out, L.rep.get(L.points[Math.floor(next() * L.points.length)]).p);
    return out;
  }

  /** Every element, one at a time. Meant for groups small enough to walk through. */
  *elements() {
    const levels = this.#levels;
    const n = this.#n;
    function* walk(i, acc) {
      if (i === levels.length) { yield acc; return; }
      for (const pt of levels[i].points) yield* walk(i + 1, mulP(acc, levels[i].rep.get(pt).p));
    }
    yield* walk(0, idP(n));
  }

  /** The orbit of a point (counted from 1), sorted. */
  orbit(point) {
    const start = point - 1;
    if (!Number.isInteger(start) || start < 0 || start >= this.#n) throw new RangeError(`${point} is not one of 1, …, ${this.#n}`);
    const seen = new Set([start]);
    const queue = [start];
    for (let h = 0; h < queue.length; h += 1) for (const g of this.#gens) { const y = g[queue[h]]; if (!seen.has(y)) { seen.add(y); queue.push(y); } }
    return queue.map((x) => x + 1).sort((a, b) => a - b);
  }

  /** All orbits on the points. */
  orbits() {
    const seen = new Set();
    const out = [];
    for (let i = 1; i <= this.#n; i += 1) if (!seen.has(i)) { const o = this.orbit(i); for (const x of o) seen.add(x); out.push(o); }
    return out;
  }

  get isTransitive() { return this.orbit(1).length === this.#n; }

  /** The subgroup generated by the given members. */
  subgroup(generators, options = {}) {
    const gens = generators.map((g) => this.#read(g));
    for (const g of gens) if (!this.contains(g)) throw new RangeError(`${cycleString(g)} is not in the group`);
    return new PermChain(gens, { degree: this.#n, ...options });
  }

  /** The subgroup fixing every one of the given points (counted from 1). */
  stabilizer(points) {
    const pts = Array.isArray(points) ? points : [points];
    const C = new PermChain(this.#gens, { degree: this.#n, base: pts, words: false });
    return new PermChain(C.#levelGenerators(pts.length), { degree: this.#n, name: null });
  }

  #levelGenerators(k) {
    const out = [];
    for (let i = k; i < this.#levels.length; i += 1) for (const g of this.#levels[i].gens) out.push(g.p);
    return out;
  }

  /** Whether H (a PermChain on the same points) lies inside this group. */
  isSubgroup(H) {
    return H.degree === this.#n && H.generators.every((g) => this.contains(g));
  }

  equals(H) { return H.degree === this.#n && H.order === this.order && this.isSubgroup(H); }

  /** The index of a subgroup, as a BigInt. */
  index(H) {
    if (!this.isSubgroup(H)) throw new RangeError("index needs a subgroup");
    return this.order / H.order;
  }

  /** Whether the subgroup H is normal here: every generator of the group carries every generator of H into H. */
  isNormal(H) {
    if (!this.isSubgroup(H)) return false;
    const hs = H.generators;
    return this.#gens.every((g) => { const gi = invP(g); return hs.every((h) => H.contains(mulP(mulP(g, h), gi))); });
  }

  /** The smallest normal subgroup holding the given members. */
  normalClosure(generators) {
    let gens = generators.map((g) => this.#read(g)).filter((g) => !isId(g));
    let N = new PermChain(gens, { degree: this.#n, words: false });
    for (;;) {
      const extra = [];
      for (const g of this.#gens) {
        const gi = invP(g);
        for (const h of gens) { const c = mulP(mulP(g, h), gi); if (!N.contains(c)) { extra.push(c); N = new PermChain([...gens, ...extra], { degree: this.#n, words: false }); } }
      }
      if (!extra.length) return new PermChain(N.strongGenerators, { degree: this.#n });
      gens = [...gens, ...extra];
    }
  }

  /** The commutator subgroup. */
  derivedSubgroup() {
    const comms = [];
    for (let i = 0; i < this.#gens.length; i += 1) for (let j = i + 1; j < this.#gens.length; j += 1) {
      const a = this.#gens[i];
      const b = this.#gens[j];
      comms.push(mulP(mulP(a, b), mulP(invP(a), invP(b))));
    }
    return this.normalClosure(comms);
  }

  /** G, G′, G″, … until it stops shrinking. */
  derivedSeries() {
    const out = [this];
    for (;;) {
      const D = out[out.length - 1].derivedSubgroup();
      if (D.order === out[out.length - 1].order) return out;
      out.push(D);
      if (D.order === 1n) return out;
    }
  }

  get isAbelian() { return this.#gens.every((a, i) => this.#gens.slice(i + 1).every((b) => samePerm(mulP(a, b), mulP(b, a)))); }
  get isSolvable() { const s = this.derivedSeries(); return s[s.length - 1].order === 1n; }
  get isPerfect() { return this.order > 0n && this.derivedSubgroup().order === this.order; }
  /** Whether the group is all of Sₙ, or all of Aₙ: decided by the order. */
  get isSymmetric() { return this.order === factorial(this.#n); }
  get isAlternating() { return this.#n >= 2 && this.order === factorial(this.#n) / 2n && this.#gens.every(evenPerm); }

  /**
   * The smallest block system in which the two points lie in one block, for a
   * transitive group: the blocks, as sorted lists of points.
   */
  minimalBlocks(a, b) {
    if (!this.isTransitive) throw new RangeError("block systems are for transitive groups");
    const n = this.#n;
    const parent = idP(n);
    const find = (x) => { while (parent[x] !== x) { parent[x] = parent[parent[x]]; x = parent[x]; } return x; };
    const queue = [];
    const merge = (x, y) => { const rx = find(x); const ry = find(y); if (rx === ry) return; const [keep, gone] = rx < ry ? [rx, ry] : [ry, rx]; parent[gone] = keep; queue.push(gone); };
    merge(a - 1, b - 1);
    while (queue.length) {
      const gamma = queue.pop();
      for (const g of this.#gens) merge(g[gamma], g[find(gamma)]);
    }
    const blocks = new Map();
    for (let i = 0; i < n; i += 1) { const r = find(i); if (!blocks.has(r)) blocks.set(r, []); blocks.get(r).push(i + 1); }
    return [...blocks.values()];
  }

  /** Every minimal block system with more than one block, each given by its blocks. */
  blockSystems() {
    if (!this.isTransitive) throw new RangeError("block systems are for transitive groups");
    const seen = new Set();
    const out = [];
    // one second point from each orbit of the stabiliser of the first is enough
    const S = this.#n > 1 ? this.stabilizer(1) : this;
    for (const o of S.orbits()) {
      if (o[0] === 1) continue;
      const blocks = this.minimalBlocks(1, o[0]);
      if (blocks.length === 1) continue;
      const id = blocks.map((b) => b.join(",")).sort().join("|");
      if (!seen.has(id)) { seen.add(id); out.push(blocks); }
    }
    return out;
  }

  /** Transitive, with no block system between single points and everything. */
  get isPrimitive() { return this.isTransitive && this.blockSystems().length === 0; }

  /** The same group with every element listed, as a `Group`; for orders up to `maxOrder`. */
  toGroup({ maxOrder = 20000, name = this.name } = {}) {
    if (this.order > BigInt(maxOrder)) throw new RangeError(`the group has ${this.order} elements, too many to list (maxOrder ${maxOrder})`);
    return permutationGroup(this.#gens.length ? this.#gens : [idP(this.#n)], this.#n, { name, maxOrder });
  }

  toString() {
    return this.name ?? `a permutation group of order ${this.order} on ${this.#n} points`;
  }

  describe() {
    return `${this.name ? `${this.name}: ` : ""}order ${this.order} = ${this.orbitLengths.join(" · ") || "1"}, base ${this.base.join(", ") || "empty"}`;
  }
}

const SUPERS = "⁰¹²³⁴⁵⁶⁷⁸⁹";
const supers = (k) => String(k).split("").map((c) => SUPERS[Number(c)]).join("");
function collapse(word) {
  const out = [];
  for (const k of word) { if (out.length && out[out.length - 1][0] === k) out[out.length - 1][1] += 1; else out.push([k, 1]); }
  return out;
}
const factorial = (n) => { let f = 1n; for (let i = 2n; i <= BigInt(n); i += 1n) f *= i; return f; };
function evenPerm(p) {
  const seen = new Uint8Array(p.length);
  let swaps = 0;
  for (let i = 0; i < p.length; i += 1) { if (seen[i]) continue; let len = 0; for (let j = i; !seen[j]; j = p[j]) { seen[j] = 1; len += 1; } swaps += len - 1; }
  return swaps % 2 === 0;
}

/** The group generated by the permutations, held by a base and strong generating set. */
export function schreierSims(generators, options = {}) {
  return new PermChain(generators, typeof options === "number" ? { degree: options } : options);
}

// ---- groups worth having at this size

const cyc = (list) => `(${list.join(" ")})`;
const rangeFrom1 = (n) => Array.from({ length: n }, (_, i) => i + 1);

/** The stickers of an N×N×N cube that turn (centres of odd cubes left out), and the six face turns as permutations of them. */
function cubeTurns(N) {
  const stickers = [];
  const c = (N - 1) / 2;
  for (let axis = 0; axis < 3; axis += 1) for (const side of [1, -1]) {
    for (let i = 0; i < N; i += 1) for (let j = 0; j < N; j += 1) {
      const pos = [0, 0, 0];
      pos[axis] = side * c;
      pos[(axis + 1) % 3] = i - c;
      pos[(axis + 2) % 3] = j - c;
      if (N % 2 === 1 && pos[(axis + 1) % 3] === 0 && pos[(axis + 2) % 3] === 0) continue;
      const normal = [0, 0, 0];
      normal[axis] = side;
      stickers.push({ pos, normal });
    }
  }
  const id = (s) => `${s.pos.join(",")}|${s.normal.join(",")}`;
  const where = new Map(stickers.map((s, i) => [id(s), i]));
  // a quarter turn about `axis` of the layer at +c or −c, clockwise seen from outside
  const turn = (axis, side) => {
    const a = (axis + 1) % 3;
    const b = (axis + 2) % 3;
    const rot = (v) => { const out = [...v]; out[a] = side * v[b]; out[b] = -side * v[a]; return out; };
    return stickers.map((s, i) => (s.pos[axis] === side * c ? where.get(id({ pos: rot(s.pos), normal: rot(s.normal) })) : i));
  };
  return { count: stickers.length, turns: { U: turn(2, 1), D: turn(2, -1), R: turn(0, 1), L: turn(0, -1), F: turn(1, -1), B: turn(1, 1) } };
}

export const largeGroups = Object.freeze({
  /** Sₙ for any n. */
  symmetric(n) {
    const gens = n < 2 ? [] : n === 2 ? ["(1 2)"] : ["(1 2)", cyc(rangeFrom1(n))];
    return new PermChain(gens, { degree: n, name: `S${n}`, names: ["s", "c"].slice(0, gens.length) });
  },
  /** Aₙ for any n. */
  alternating(n) {
    const gens = n < 3 ? [] : n % 2 ? ["(1 2 3)", cyc(rangeFrom1(n))] : ["(1 2 3)", cyc(rangeFrom1(n).slice(1))];
    return new PermChain(gens, { degree: Math.max(n, 1), name: `A${n}` });
  },
  /** The Mathieu groups M₁₁, M₁₂, M₂₂, M₂₃ and M₂₄. */
  mathieu(n) {
    const a23 = cyc(rangeFrom1(23));
    const b23 = "(3 17 10 7 9)(4 13 14 19 5)(8 18 11 12 23)(15 20 22 21 16)";
    const c24 = "(1 24)(2 23)(3 12)(4 16)(5 18)(6 10)(7 20)(8 14)(9 21)(11 17)(13 22)(15 19)";
    if (n === 11) return new PermChain([cyc(rangeFrom1(11)), "(3 7 11 8)(4 10 5 6)"], { degree: 11, name: "M11" });
    if (n === 12) return new PermChain([cyc(rangeFrom1(11)), "(3 7 11 8)(4 10 5 6)", "(1 12)(2 11)(3 6)(4 8)(5 9)(7 10)"], { degree: 12, name: "M12" });
    if (n === 23) return new PermChain([a23, b23], { degree: 23, name: "M23" });
    if (n === 24) return new PermChain([a23, b23, c24], { degree: 24, name: "M24" });
    if (n === 22) {
      const S = new PermChain([a23, b23], { degree: 23 }).stabilizer(23);
      return new PermChain(S.generators.map((g) => g.slice(0, 22)), { degree: 22, name: "M22" });
    }
    throw new RangeError("the Mathieu groups are M11, M12, M22, M23 and M24");
  },
  /** The group of Rubik's cube: the six face turns acting on the 48 stickers that move. */
  rubik() {
    const { count, turns } = cubeTurns(3);
    const names = ["U", "D", "L", "R", "F", "B"];
    return new PermChain(names.map((k) => turns[k]), { degree: count, names, name: "Rubik's cube group" });
  },
  /** The 2×2×2 cube with one corner held still: the turns U, R, F on its 24 stickers. */
  pocketCube() {
    const { count, turns } = cubeTurns(2);
    const names = ["U", "R", "F"];
    return new PermChain(names.map((k) => turns[k]), { degree: count, names, name: "pocket cube group" });
  },
});
