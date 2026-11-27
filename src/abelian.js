// Finitely generated abelian groups (and Z-modules), through the Smith
// normal form.
//
// Generators x₁ … xₙ with relations r₁, …, r_m (each a row of integers,
// meaning Σ r_ij x_j = 0) give the group Zⁿ/L, L the row lattice. If
// U·R·V = D is the Smith form, the new generators y = V⁻¹x satisfy
// dᵢ yᵢ = 0 and nothing else, so the group is Z/d₁ × … × Z/d_k × Z^(n−k)
// with the units dᵢ = 1 dropped: the structure theorem, computed.

import { smithForm, ZZ } from "./smith.js";
import { factorInt } from "./group.js";
import { cyclic, directProduct } from "./named.js";

const big = (x) => ZZ.from(x);
const fmt = (d) => (d === 0n ? "Z" : `Z/${d}`);

export class AbelianGroup {
  /**
   * The group with the given invariant factors (d₁ | d₂ | …, each > 1) and
   * free rank. `generators` (optional) says how each cyclic factor is made
   * from the original generators, as integer vectors.
   */
  constructor(invariants, rank = 0, { generators = null, names = null } = {}) {
    const inv = invariants.map(big);
    for (let i = 0; i < inv.length; i += 1) {
      if (inv[i] <= 1n) throw new RangeError("invariant factors are whole numbers above 1");
      if (i && inv[i] % inv[i - 1] !== 0n) throw new RangeError(`invariant factors must divide each other: ${inv[i - 1]} does not divide ${inv[i]}`);
    }
    this.invariants = Object.freeze(inv);
    this.rank = rank;
    this.generators = generators;
    this.names = names;
    Object.freeze(this);
  }

  /** The order, or null when the group is infinite. */
  get order() { return this.rank ? null : this.invariants.reduce((a, b) => a * b, 1n); }
  get isFinite() { return this.rank === 0; }
  get isTrivial() { return this.rank === 0 && this.invariants.length === 0; }
  get isCyclic() { return this.invariants.length + this.rank <= 1; }
  /** The torsion subgroup's invariants and order. */
  get torsion() { return new AbelianGroup(this.invariants, 0); }
  /** The exponent of the torsion part (the largest invariant factor), 1 if none. */
  get exponent() { return this.invariants.length ? this.invariants[this.invariants.length - 1] : 1n; }

  /** The prime-power cyclic factors: [{ p, power, order }], grouped by p. */
  get elementaryDivisors() {
    const out = [];
    for (const d of this.invariants) {
      for (const [p, k] of factorInt(Number(d))) out.push({ p: BigInt(p), power: k, order: BigInt(p) ** BigInt(k) });
    }
    return out.sort((a, b) => (a.p === b.p ? a.power - b.power : a.p < b.p ? -1 : 1));
  }

  /** Same group? (Same invariants and rank: the structure theorem says that is enough.) */
  equals(H) { return H.rank === this.rank && H.invariants.length === this.invariants.length && H.invariants.every((d, i) => d === this.invariants[i]); }

  /** G × H. */
  times(H) {
    const ed = [...this.elementaryDivisors, ...H.elementaryDivisors].map((e) => e.order);
    return fromElementaryDivisors(ed, this.rank + H.rank);
  }

  /** The finite group as a Group: a product of cyclic groups (order at most 20,000). */
  toGroup() {
    if (!this.isFinite) throw new RangeError("an infinite group cannot be listed");
    if (!this.invariants.length) return cyclic(1);
    const cs = this.invariants.map((d) => cyclic(Number(d)));
    return cs.length === 1 ? cs[0] : directProduct(...cs, { name: this.toString() });
  }

  /** "Z/2 × Z/12 × Z²", or "1" for the trivial group. With { primary: true }, by elementary divisors. */
  toString({ primary = false } = {}) {
    const parts = primary ? this.elementaryDivisors.map((e) => fmt(e.order)) : this.invariants.map(fmt);
    if (this.rank) parts.push(this.rank === 1 ? "Z" : `Z^${this.rank}`);
    return parts.length ? parts.join(" × ") : "1";
  }
}

/** The abelian group with these invariant factors and free rank (units are dropped, the list is checked). */
export const abelianGroup = (invariants = [], rank = 0) => new AbelianGroup(invariants.map(big).filter((d) => d !== 1n), rank);

/** Z/n₁ × Z/n₂ × … × Z^rank for any orders nᵢ (0 means Z), put into invariant factor form. */
export function cyclicProduct(orders, rank = 0) {
  const o = orders.map(big);
  if (o.some((d) => d < 0n)) throw new RangeError("orders are whole numbers ≥ 0 (0 for Z)");
  const rows = [];
  o.forEach((d, i) => { if (d !== 0n) rows.push(o.map((_, j) => (i === j ? d : 0n))); });
  return fromRelations(rows, o.length, rank);
}

function fromElementaryDivisors(orders, rank) {
  // group the prime powers by prime, largest last, and multiply across
  const byP = new Map();
  for (const q of orders) {
    for (const [p] of factorInt(Number(q))) {
      if (!byP.has(p)) byP.set(p, []);
      byP.get(p).push(q);
    }
  }
  let k = 0;
  for (const l of byP.values()) { l.sort((a, b) => (a < b ? -1 : 1)); k = Math.max(k, l.length); }
  const inv = Array(k).fill(1n);
  for (const l of byP.values()) l.forEach((q, i) => { inv[k - l.length + i] *= q; });
  return new AbelianGroup(inv.filter((d) => d !== 1n), rank);
}

/**
 * The abelian group with generators x₁ … xₙ and the given relations (rows of
 * integers: Σ r_j x_j = 0). Unlisted relations mean free directions.
 * The result also says, for each cyclic factor, which combination of the
 * original generators generates it (`generators`, integer vectors).
 */
export function fromRelations(relations, n = relations[0]?.length ?? 0, extraRank = 0) {
  const rows = relations.length ? relations.map((r) => {
    if (r.length !== n) throw new RangeError(`each relation needs ${n} entries`);
    return r.map(big);
  }) : [Array(n).fill(0n)];
  const S = smithForm(rows, ZZ);
  const d = [...S.diagonal];
  while (d.length < n) d.push(0n);
  const invariants = [];
  const gens = [];
  let rank = extraRank;
  d.forEach((di, i) => {
    if (di === 1n) return;
    if (di === 0n) rank += 1;
    else invariants.push(di);
    gens.push({ order: di === 0n ? null : di, vector: S.Vinv[i] });
  });
  // free generators after the torsion ones
  gens.sort((a, b) => (a.order === null) - (b.order === null));
  return new AbelianGroup(invariants, rank, { generators: gens });
}

/** Every abelian group of order n, one of each kind (as many as Π p(aᵢ) for n = Π pᵢ^aᵢ). */
export function abelianGroups(n) {
  if (!Number.isInteger(n) || n < 1) throw new RangeError("n must be a positive whole number");
  const parts = (a, max = a) => {
    if (a === 0) return [[]];
    const out = [];
    for (let k = Math.min(a, max); k >= 1; k -= 1) for (const rest of parts(a - k, k)) out.push([k, ...rest]);
    return out;
  };
  let choices = [[]];
  for (const [p, a] of factorInt(n)) {
    const next = [];
    for (const c of choices) for (const lam of parts(a)) next.push([...c, ...lam.map((k) => BigInt(p) ** BigInt(k))]);
    choices = next;
  }
  return choices.map((orders) => fromElementaryDivisors(orders, 0));
}
