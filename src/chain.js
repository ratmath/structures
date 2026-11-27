// Chain complexes of finitely generated free abelian groups, and their
// homology.
//
// A chain complex is a sequence of free groups C_k = Z^(n_k) (one generator
// per k-cell) and boundary maps ∂_k : C_k → C_(k−1) with ∂_(k−1)∘∂_k = 0.
// Homology H_k = ker ∂_k / im ∂_(k+1) measures the cycles that are not
// boundaries. Everything comes from the Smith normal forms of the ∂_k: if
// ∂_k has rank r_k and nonzero invariant factors d_1 | … | d_(r_k), then
//
//   H_k ≅ Z^(n_k − r_k − r_(k+1)) × Z/d_i for each d_i > 1 of ∂_(k+1),
//
// over Q the dimension is n_k − r_k − r_(k+1), and over Z/p the ranks are
// counted mod p: the number of invariant factors not divisible by p (the
// Smith form is reached by invertible integer operations, which stay
// invertible mod p).
//
// Matrices follow the usual convention: ∂_k has one column per k-cell and one
// row per (k−1)-cell, the column of a cell holding its boundary.

import { smithForm, ZZ } from "./smith.js";
import { abelianGroup } from "./abelian.js";

const big = (x) => ZZ.from(x);

const isPrime = (p) => {
  if (p < 2n) return false;
  for (let d = 2n; d * d <= p; d += 1n) if (p % d === 0n) return false;
  return true;
};

/** Homology over a field: a vector space known by its dimension. */
class FieldHomology {
  constructor(field, dimension) {
    this.field = field;
    this.dimension = dimension;
    Object.freeze(this);
  }
  get rank() { return this.dimension; }
  get isTrivial() { return this.dimension === 0; }
  toString() {
    if (!this.dimension) return "0";
    const F = this.field === "Q" ? "Q" : `(${this.field})`;
    return this.dimension === 1 ? this.field : `${F}^${this.dimension}`;
  }
}

export class ChainComplex {
  /**
   * boundaries[k − 1] is ∂_k : C_k → C_(k−1), so the list starts with ∂_1.
   * sizes (optional) gives n_0, n_1, … when a matrix has no rows or columns.
   */
  constructor(boundaries, { sizes = null, check = true } = {}) {
    if (!Array.isArray(boundaries)) throw new TypeError("chainComplex needs a list of boundary matrices [∂1, ∂2, ...]");
    const mats = boundaries.map((M) => M.map((row) => row.map(big)));
    const top = Math.max(mats.length, sizes ? sizes.length - 1 : 0);
    const cols = (M) => (M && M.length ? M[0].length : null);
    const n = [];
    for (let k = 0; k <= top; k += 1) {
      const given = sizes ? sizes[k] ?? null : null;
      const below = k >= 1 ? cols(mats[k - 1]) : null; // columns of ∂_k
      const above = mats[k] && mats[k].length ? mats[k].length : null; // rows of ∂_(k+1)
      const nk = given ?? below ?? above ?? (k === 0 && mats[0] ? 0 : null);
      if (nk === null) throw new RangeError(`the number of ${k}-cells cannot be read from empty matrices: pass { sizes }`);
      for (const v of [below, above]) if (v !== null && v !== nk) throw new RangeError(`C_${k} has ${nk} cells from one matrix and ${v} from another`);
      n.push(nk);
    }
    for (let k = 1; k <= mats.length; k += 1) {
      const M = mats[k - 1];
      if (M.some((r) => r.length !== M[0].length)) throw new RangeError(`the rows of ∂${k} have different lengths`);
      if (M.length && M.length !== n[k - 1]) throw new RangeError(`∂${k} needs ${n[k - 1]} rows (one per ${k - 1}-cell), not ${M.length}`);
      if (!M.length && n[k - 1] && n[k]) throw new RangeError(`∂${k} is empty but there are ${n[k - 1]} ${k - 1}-cells and ${n[k]} ${k}-cells`);
    }
    this.sizes = Object.freeze(n);
    this.dimension = n.length - 1;
    // pad to full shape: ∂_k is n_(k−1) × n_k
    this.boundaries = Object.freeze(n.slice(1).map((nk, i) => {
      const M = mats[i];
      if (M && M.length && M[0].length) return Object.freeze(M.map((r) => Object.freeze(r)));
      return Object.freeze(Array.from({ length: n[i] }, () => Object.freeze(Array(nk).fill(0n))));
    }));
    if (check) {
      for (let k = 2; k <= this.dimension; k += 1) {
        const A = this.boundaries[k - 2]; // ∂_(k−1): n_(k−2) × n_(k−1)
        const B = this.boundaries[k - 1]; // ∂_k: n_(k−1) × n_k
        for (let i = 0; i < A.length; i += 1) for (let j = 0; j < n[k]; j += 1) {
          let s = 0n;
          for (let t = 0; t < n[k - 1]; t += 1) s += A[i][t] * B[t][j];
          if (s !== 0n) throw new RangeError(`not a chain complex: ∂${k - 1}∘∂${k} is not zero (entry ${i + 1}, ${j + 1} is ${s})`);
        }
      }
    }
    this._smith = new Map();
    Object.freeze(this);
  }

  /** ∂_k as a matrix (n_(k−1) rows, n_k columns); the zero map outside 1 … dimension. */
  boundary(k) {
    if (k >= 1 && k <= this.dimension) return this.boundaries[k - 1];
    return null;
  }

  /** The nonzero invariant factors of ∂_k (empty for the zero map). */
  invariantFactors(k) {
    if (k < 1 || k > this.dimension) return [];
    if (!this._smith.has(k)) {
      const M = this.boundaries[k - 1];
      const d = M.length && M[0].length ? smithForm(M, ZZ).diagonal.filter((x) => x !== 0n) : [];
      this._smith.set(k, Object.freeze(d));
    }
    return this._smith.get(k);
  }

  /** The rank of ∂_k over Z (= over Q), or over Z/p when p is given. */
  rank(k, p = null) {
    const d = this.invariantFactors(k);
    if (p === null) return d.length;
    return d.filter((x) => x % p !== 0n).length;
  }

  /** Normalises a coefficient choice to null (Z), "Q", or a prime BigInt. */
  _field(coefficients) {
    if (coefficients === undefined || coefficients === null || coefficients === "Z") return null;
    if (coefficients === "Q") return "Q";
    const m = typeof coefficients === "string" ? coefficients.match(/^Z\/(\d+)$/) : null;
    const p = m ? BigInt(m[1]) : typeof coefficients === "number" || typeof coefficients === "bigint" ? BigInt(coefficients) : null;
    if (p === null || !isPrime(p)) throw new RangeError(`coefficients are "Z", "Q" or "Z/p" for a prime p, not ${coefficients}`);
    return p;
  }

  /**
   * H_k. Over Z (the default) an AbelianGroup such as "Z^2" or "Z × Z/2";
   * with { coefficients: "Q" } or "Z/p" a vector space known by its dimension.
   */
  homology(k, { coefficients = "Z" } = {}) {
    if (!Number.isInteger(k) || k < 0) throw new RangeError("homology(k) needs a whole number k ≥ 0");
    if (k > this.dimension) return coefficients === "Z" || !coefficients ? abelianGroup([], 0) : new FieldHomology(String(coefficients), 0);
    const F = this._field(coefficients);
    const nk = this.sizes[k];
    if (F === null) {
      const free = nk - this.rank(k) - this.rank(k + 1);
      const torsion = this.invariantFactors(k + 1).filter((d) => d > 1n);
      return abelianGroup(torsion, free);
    }
    const p = F === "Q" ? null : F;
    return new FieldHomology(F === "Q" ? "Q" : `Z/${F}`, nk - this.rank(k, p) - this.rank(k + 1, p));
  }

  /** Every H_k for k = 0 … dimension. */
  homologyGroups(options = {}) {
    return Array.from({ length: this.dimension + 1 }, (_, k) => this.homology(k, options));
  }

  /** The Betti numbers b_k = rank H_k (over Q, or over Z/p with { coefficients }). */
  betti(options = {}) {
    return this.homologyGroups(options).map((H) => H.rank);
  }

  /** The torsion of each H_k over Z, as AbelianGroups. */
  get torsion() {
    return this.homologyGroups().map((H) => H.torsion);
  }

  /** χ = Σ (−1)^k n_k, which equals Σ (−1)^k b_k over any field. */
  get eulerCharacteristic() {
    return this.sizes.reduce((s, nk, k) => s + (k % 2 ? -nk : nk), 0);
  }

  /** The reduced homology: H̃_0 has rank one less (for a nonempty complex); the rest are unchanged. */
  reducedHomology(k, options = {}) {
    const H = this.homology(k, options);
    if (k !== 0 || this.sizes[0] === 0) return H;
    if (H instanceof FieldHomology) return new FieldHomology(H.field, H.dimension - 1);
    return abelianGroup(H.invariants, H.rank - 1);
  }

  /** "Z, Z^2, Z" style summary of the integral homology. */
  toString() {
    return this.homologyGroups().map((H, k) => `H${k} = ${H instanceof FieldHomology ? H : H.toString({ additive: true })}`).join(", ");
  }
}

/**
 * A chain complex from its boundary matrices [∂1, ∂2, …], ∂_k with a column
 * per k-cell and a row per (k−1)-cell. Checks ∂∂ = 0 and refuses otherwise.
 * { sizes: [n0, n1, …] } fixes the number of cells where a map is empty.
 *
 *   chainComplex([[[0, 0]], [[0], [0]]])        // the torus: H = Z, Z^2, Z
 *   chainComplex([[[0]], [[2]]])                // RP²: H = Z, Z/2, 0
 */
export function chainComplex(boundaries, options = {}) {
  return new ChainComplex(boundaries, options);
}
