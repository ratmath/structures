// Group actions on finite sets: orbits, stabilisers, fixed points and the
// orbit-counting lemma.
//
// An action of G on a set X is act(g, x) with act(e, x) = x and
// act(gh, x) = act(g, act(h, x)). Each group element becomes a permutation of
// X; those permutations are built from the generators along a spanning tree
// of G and checked against act, so a function that is not an action is
// caught rather than giving wrong counts.

import { PermGroup } from "@ratmath/algebraic";
import { Rational } from "@ratmath/core";
import { defaultKey } from "./group.js";

const CHECK_LIMIT = 2_000_000;
const POINT_LIMIT = 1_000_000;

export class Action {
  #perms;
  #index;
  #key;

  constructor(group, points, act, { key = defaultKey, format = null, name = null, check } = {}) {
    if (points.length > POINT_LIMIT) throw new RangeError(`an action on more than ${POINT_LIMIT} points`);
    this.group = group;
    this.points = Object.freeze([...points]);
    this.act = act;
    this.name = name;
    this.#key = key;
    this.formatPoint = format ?? ((x) => (typeof x === "string" ? x : Array.isArray(x) ? `(${x.join(",")})` : String(x)));
    this.#index = new Map(this.points.map((x, i) => [key(x), i]));
    if (this.#index.size !== this.points.length) throw new RangeError("two points have the same key");
    this.#build(check ?? group.order * points.length <= CHECK_LIMIT);
  }

  get degree() { return this.points.length; }

  pointIndex(x) {
    const i = this.#index.get(this.#key(x));
    if (i === undefined) throw new RangeError(`${this.formatPoint(x)} is not one of the points`);
    return i;
  }

  #build(check) {
    const G = this.group;
    const m = this.points.length;
    const genPerm = [...G.gens].map((s) => {
      const g = G.elements[s];
      const p = new Int32Array(m);
      this.points.forEach((x, i) => {
        const y = this.#index.get(this.#key(this.act(g, x)));
        if (y === undefined) throw new RangeError(`${G.format(g)} sends ${this.formatPoint(x)} outside the set of points`);
        p[i] = y;
      });
      if (new Set(p).size !== m) throw new RangeError(`${G.format(g)} does not act one-to-one, so this is not an action`);
      return p;
    });
    // x = parent·s acts as perm(parent)∘perm(s): first s, then parent
    const perms = new Array(G.order);
    perms[0] = Int32Array.from({ length: m }, (_, i) => i);
    const { order, from, by } = G.spanningTree();
    for (const x of order) {
      if (x === 0) continue;
      const a = perms[from[x]];
      const s = genPerm[by[x]];
      const p = new Int32Array(m);
      for (let i = 0; i < m; i += 1) p[i] = a[s[i]];
      perms[x] = p;
    }
    // the same rule must hold along every edge, or act is not compatible with the group's multiplication
    for (let x = 0; x < G.order; x += 1) {
      for (let k = 0; k < G.gens.length; k += 1) {
        const y = G.mulIndex(x, G.gens[k]);
        const a = perms[x];
        const s = genPerm[k];
        const p = perms[y];
        for (let i = 0; i < m; i += 1) if (p[i] !== a[s[i]]) throw new RangeError("act(gh, x) ≠ act(g, act(h, x)) somewhere, so this is not an action");
      }
    }
    if (check) {
      for (let g = 0; g < G.order; g += 1) {
        const p = perms[g];
        for (let i = 0; i < m; i += 1) {
          if (this.#index.get(this.#key(this.act(G.elements[g], this.points[i]))) !== p[i]) throw new RangeError(`act disagrees with the group's multiplication at ${G.format(G.elements[g])}, so this is not an action`);
        }
      }
    }
    this.#perms = perms;
  }

  /** The permutation of point numbers by which g acts. */
  permutation(g) { return [...this.#perms[this.group.indexOf(g) < 0 ? this.#bad(g) : this.group.indexOf(g)]]; }
  #bad(g) { throw new RangeError(`${g} is not in the group`); }
  permutationIndex(i) { return this.#perms[i]; }

  /** g·x. */
  apply(g, x) { return this.points[this.#perms[this.group.indexOf(g)][this.pointIndex(x)]]; }

  /** The orbit G·x, in the order it is reached. */
  orbit(x) {
    const start = this.pointIndex(x);
    const seen = new Uint8Array(this.points.length);
    seen[start] = 1;
    const list = [start];
    const gp = [...this.group.gens].map((s) => this.#perms[s]);
    for (let i = 0; i < list.length; i += 1) for (const p of gp) {
      const y = p[list[i]];
      if (!seen[y]) { seen[y] = 1; list.push(y); }
    }
    return list.map((i) => this.points[i]);
  }

  /** All orbits, each a list of points; together they partition the points. */
  orbits() {
    const seen = new Uint8Array(this.points.length);
    const out = [];
    const gp = [...this.group.gens].map((s) => this.#perms[s]);
    for (let s = 0; s < this.points.length; s += 1) {
      if (seen[s]) continue;
      seen[s] = 1;
      const list = [s];
      for (let i = 0; i < list.length; i += 1) for (const p of gp) {
        const y = p[list[i]];
        if (!seen[y]) { seen[y] = 1; list.push(y); }
      }
      out.push(list.map((i) => this.points[i]));
    }
    return out;
  }

  /** Stab(x) = { g : g·x = x }, a subgroup. |orbit|·|Stab| = |G|. */
  stabiliser(x) {
    const i = this.pointIndex(x);
    const list = [];
    for (let g = 0; g < this.group.order; g += 1) if (this.#perms[g][i] === i) list.push(g);
    return this.group.subgroupFromIndices(list);
  }
  stabilizer(x) { return this.stabiliser(x); }

  /** Fix(g) = { x : g·x = x }. */
  fixedPoints(g) {
    const p = this.#perms[this.group.indexOf(g)];
    if (!p) throw new RangeError(`${g} is not in the group`);
    return this.points.filter((_, i) => p[i] === i);
  }

  /** The kernel: the elements that fix every point (a normal subgroup). */
  kernel() {
    const list = [];
    for (let g = 0; g < this.group.order; g += 1) if (this.#perms[g].every((y, i) => y === i)) list.push(g);
    return this.group.subgroupFromIndices(list);
  }

  get isTransitive() { return this.points.length <= 1 || this.orbit(this.points[0]).length === this.points.length; }
  get isFaithful() { return this.kernel().order === 1; }
  /** Free: only the identity fixes any point. */
  get isFree() {
    for (let g = 1; g < this.group.order; g += 1) if (this.#perms[g].some((y, i) => y === i)) return false;
    return true;
  }

  /**
   * The orbit-counting lemma (Burnside, Cauchy–Frobenius):
   * number of orbits = (1/|G|) Σ_g |Fix(g)|. Returns
   * { count, total, order, byClass } with total = Σ |Fix(g)| and byClass
   * the fixed-point count for each conjugacy class (constant on a class).
   * The count is checked against the orbits found directly.
   */
  countOrbits() {
    const G = this.group;
    let total = 0n;
    const fix = new Int32Array(G.order);
    for (let g = 0; g < G.order; g += 1) {
      const p = this.#perms[g];
      let f = 0;
      for (let i = 0; i < p.length; i += 1) if (p[i] === i) f += 1;
      fix[g] = f;
      total += BigInt(f);
    }
    const q = new Rational(total, BigInt(G.order));
    if (q.denominator !== 1n) throw new Error("internal: Burnside's sum is not a multiple of |G|");
    const count = Number(q.numerator);
    const direct = this.orbits().length;
    if (count !== direct) throw new Error(`internal: Burnside gives ${count} orbits but ${direct} were found`);
    const byClass = G.conjugacyClasses().map((c) => ({ representative: c.representative, size: c.size, fixed: fix[c.indices[0]] }));
    return Object.freeze({ count, total, order: G.order, byClass, toString: () => `(1/${G.order})·(${byClass.map((c) => (c.size === 1 ? `${c.fixed}` : `${c.size}·${c.fixed}`)).join(" + ")}) = ${count}` });
  }

  /** The image of G in the symmetric group on the points, as algebraic's PermGroup on 0 … m−1 (for discrete's cycleIndex). */
  toPermGroup() { return new PermGroup([...this.group.gens].map((s) => [...this.#perms[s]]), this.points.length); }

  /** The induced action on the k-element subsets of the points (as sorted lists of point numbers' points). */
  onSubsets(k) {
    const m = this.points.length;
    const subsets = [];
    const rec = (s, a) => {
      if (a.length === k) { subsets.push(a.slice()); return; }
      for (let i = s; i < m; i += 1) { a.push(i); rec(i + 1, a); a.pop(); }
    };
    rec(0, []);
    const pts = subsets.map((s) => s.map((i) => this.points[i]));
    const keyOf = (S) => S.map((x) => this.pointIndex(x)).sort((a, b) => a - b).join(",");
    return new Action(this.group, pts, (g, S) => S.map((x) => this.apply(g, x)), { key: (S) => (Array.isArray(S) ? keyOf(S) : String(S)), format: (S) => `{${S.map((x) => this.formatPoint(x)).join(", ")}}`, name: `${this.name ?? "action"} on ${k}-subsets`, check: false });
  }

  /**
   * The induced action on colourings: every way of giving each point one of
   * `colours` (a count k, or a list of names), with (g·c)(g·x) = c(x). Its
   * orbits are the colourings up to symmetry.
   */
  onColourings(colours) {
    const names = typeof colours === "number" ? [...Array(colours).keys()] : [...colours];
    const k = names.length;
    const m = this.points.length;
    if (k ** m > POINT_LIMIT) throw new RangeError(`${k}^${m} colourings is too many to list; count them with discrete's countColourings(action.toPermGroup(), k) instead`);
    const pts = [];
    for (let v = 0; v < k ** m; v += 1) {
      const c = new Array(m);
      let t = v;
      for (let i = 0; i < m; i += 1) { c[i] = names[t % k]; t = Math.floor(t / k); }
      pts.push(c);
    }
    const self = this;
    return new Action(this.group, pts, (g, c) => {
      const p = self.#perms[self.group.indexOf(g)];
      const out = new Array(m);
      for (let i = 0; i < m; i += 1) out[p[i]] = c[i];
      return out;
    }, { key: (c) => c.join("\u0001"), format: (c) => c.join(""), name: `${this.name ?? "action"} on ${k}-colourings`, check: false });
  }
  onColorings(colours) { return this.onColourings(colours); }

  toString() { return this.name ?? `an action of ${this.group} on ${this.points.length} points`; }
}

/** An action of G on a list of points by act(g, x). */
export const action = (group, points, act, options) => new Action(group, points, act, options);

/** G acting on itself by left multiplication (Cayley's theorem). */
export const regularAction = (G) => new Action(G, G.elements, (g, x) => G.mul(g, x), { key: (x) => G.key(x), format: (x) => G.format(x), name: `${G} acting on itself`, check: false });

/** G acting on itself by conjugation; the orbits are the conjugacy classes. */
export const conjugationAction = (G) => new Action(G, G.elements, (g, x) => G.conjugate(x, g), { key: (x) => G.key(x), format: (x) => G.format(x), name: `${G} by conjugation`, check: false });

/** G acting on the left cosets of H by left multiplication; the points are coset representatives' cosets. */
export function cosetAction(G, H) {
  const cs = G.cosets(H);
  const of = new Int32Array(G.order);
  cs.forEach((c, i) => c.indices.forEach((x) => { of[x] = i; }));
  const pts = cs.map((_, i) => i);
  const act = (g, i) => of[G.mulIndex(G.indexOf(g), cs[i].indices[0])];
  const A = new Action(G, pts, act, { format: (i) => `${G.format(cs[i].representative)}H`, name: `${G} on the cosets of a subgroup of order ${H.order}`, check: false });
  A.cosets = cs;
  return A;
}

/** The orbit-counting lemma for an action: see Action#countOrbits. */
export const burnside = (A) => A.countOrbits();
