// The Smith normal form over a Euclidean domain: the integers or the
// polynomials in one variable over a field.
//
// For a matrix A it finds invertible U and V (over the domain itself, so with
// determinant a unit) with U·A·V = D diagonal and d₁ | d₂ | … . The dᵢ, made
// canonical (positive integers, monic polynomials), are the invariant
// factors; they decide the structure of a finitely generated module over the
// domain. V⁻¹ is tracked as well, because it says how the new basis is made
// from the old one.

import { ring } from "./mpoly.js";

const absB = (a) => (a < 0n ? -a : a);

/** The integers as a Euclidean domain (elements are BigInts; numbers are accepted). */
export const ZZ = Object.freeze({
  name: "Z",
  zero: 0n,
  one: 1n,
  from: (x) => {
    if (typeof x === "bigint") return x;
    if (typeof x === "number" && Number.isInteger(x)) return BigInt(x);
    if (typeof x === "string" && /^\s*-?\d+\s*$/.test(x)) return BigInt(x.trim());
    if (x && typeof x.numerator === "bigint" && x.denominator === 1n) return x.numerator;
    throw new TypeError(`${x} is not a whole number`);
  },
  isZero: (a) => a === 0n,
  add: (a, b) => a + b,
  sub: (a, b) => a - b,
  mul: (a, b) => a * b,
  neg: (a) => -a,
  /** Euclidean division with 0 ≤ r < |b|. */
  divmod: (a, b) => {
    let q = a / b;
    let r = a - q * b;
    if (r < 0n) { r += absB(b); q = (a - r) / b; }
    return { quotient: q, remainder: r };
  },
  norm: (a) => absB(a),
  /** The unit u with u·a canonical (positive), and its inverse. */
  unit: (a) => (a < 0n ? { u: -1n, inv: -1n } : { u: 1n, inv: 1n }),
  equals: (a, b) => a === b,
  format: (a) => String(a),
});

/**
 * The polynomials in one variable over Q (or over a field object) as a
 * Euclidean domain; elements are structures polynomials. Q[x] by default.
 */
export function polynomialDomain(v = "x", options = {}) {
  const R = ring(v, options);
  return Object.freeze({
    name: `${R.field?.name ?? "Q"}[${R.vars[0]}]`,
    ring: R,
    zero: R.zero,
    one: R.one,
    from: (x) => R.from(x),
    isZero: (a) => a.isZero(),
    add: (a, b) => a.add(b),
    sub: (a, b) => a.sub(b),
    mul: (a, b) => a.mul(b),
    neg: (a) => a.neg(),
    divmod: (a, b) => a.divmod(b),
    norm: (a) => (a.isZero() ? -1 : a.degree()),
    unit: (a) => {
      if (a.isZero()) return { u: R.one, inv: R.one };
      const c = R._coeff(a.leadingCoefficient());
      return { u: R.one._scale(R.K.inv(c)), inv: R.one._scale(c) };
    },
    equals: (a, b) => a.equals(b),
    format: (a) => String(a),
  });
}

const less = (D, a, b) => {
  const x = D.norm(a);
  const y = D.norm(b);
  return x < y;
};

/**
 * Smith normal form of a matrix (an array of rows) over the domain D
 * (default the integers). Returns { D, U, V, Vinv, diagonal, rank } with
 * U·A·V = D, V·Vinv = I, and diagonal the canonical d₁ | d₂ | … (as many
 * as min(rows, columns); zeros at the end).
 */
export function smithForm(A0, domain = ZZ) {
  const R = domain;
  const m = A0.length;
  const n = m ? A0[0].length : 0;
  if (A0.some((r) => r.length !== n)) throw new RangeError("the rows have different lengths");
  const A = A0.map((r) => r.map((x) => R.from(x)));
  const id = (k) => Array.from({ length: k }, (_, i) => Array.from({ length: k }, (_, j) => (i === j ? R.one : R.zero)));
  const U = id(m);
  const V = id(n);
  const Vinv = id(n);

  // elementary operations, mirrored on U (rows), V (columns) and V⁻¹ (rows)
  const swapRows = (i, j) => { [A[i], A[j]] = [A[j], A[i]]; [U[i], U[j]] = [U[j], U[i]]; };
  const addRow = (i, j, q) => { // row_i += q·row_j
    for (let k = 0; k < n; k += 1) A[i][k] = R.add(A[i][k], R.mul(q, A[j][k]));
    for (let k = 0; k < m; k += 1) U[i][k] = R.add(U[i][k], R.mul(q, U[j][k]));
  };
  const scaleRow = (i, u) => {
    for (let k = 0; k < n; k += 1) A[i][k] = R.mul(u, A[i][k]);
    for (let k = 0; k < m; k += 1) U[i][k] = R.mul(u, U[i][k]);
  };
  const swapCols = (i, j) => {
    for (const M of [A, V]) for (const row of M) [row[i], row[j]] = [row[j], row[i]];
    [Vinv[i], Vinv[j]] = [Vinv[j], Vinv[i]];
  };
  const addCol = (i, j, q) => { // col_i += q·col_j; then V⁻¹: row_j −= q·row_i
    for (const M of [A, V]) for (const row of M) row[i] = R.add(row[i], R.mul(q, row[j]));
    for (let k = 0; k < n; k += 1) Vinv[j][k] = R.sub(Vinv[j][k], R.mul(q, Vinv[i][k]));
  };

  const t0 = Math.min(m, n);
  let t = 0;
  for (; t < t0; t += 1) {
    // the nonzero entry of least norm in the lower right block becomes the pivot
    for (;;) {
      let best = null;
      for (let i = t; i < m; i += 1) for (let j = t; j < n; j += 1) {
        if (R.isZero(A[i][j])) continue;
        if (!best || less(R, A[i][j], A[best[0]][best[1]])) best = [i, j];
      }
      if (!best) break;
      if (best[0] !== t) swapRows(t, best[0]);
      if (best[1] !== t) swapCols(t, best[1]);
      const p = A[t][t];
      let clean = true;
      for (let i = t + 1; i < m; i += 1) {
        if (R.isZero(A[i][t])) continue;
        const { quotient } = R.divmod(A[i][t], p);
        addRow(i, t, R.neg(quotient));
        if (!R.isZero(A[i][t])) clean = false;
      }
      for (let j = t + 1; j < n; j += 1) {
        if (R.isZero(A[t][j])) continue;
        const { quotient } = R.divmod(A[t][j], p);
        addCol(j, t, R.neg(quotient));
        if (!R.isZero(A[t][j])) clean = false;
      }
      if (!clean) continue;
      // the pivot must divide everything left; if not, bring an offending row up and go again
      let bad = -1;
      for (let i = t + 1; i < m && bad < 0; i += 1) for (let j = t + 1; j < n; j += 1) {
        if (!R.isZero(R.divmod(A[i][j], p).remainder)) { bad = i; break; }
      }
      if (bad < 0) break;
      addRow(t, bad, R.one);
    }
    if (R.isZero(A[t][t])) break;
    const { u } = R.unit(A[t][t]);
    if (!R.equals(u, R.one)) scaleRow(t, u);
  }
  const diagonal = [];
  for (let i = 0; i < t0; i += 1) diagonal.push(A[i][i]);
  const rank = diagonal.filter((d) => !R.isZero(d)).length;
  return Object.freeze({ D: A, U, V, Vinv, diagonal, rank, domain: R });
}

/** Matrix product over a domain (used by the checks and the tests). */
export function matMul(R, X, Y) {
  return X.map((row) => Y[0].map((_, j) => row.reduce((s, x, k) => R.add(s, R.mul(x, Y[k][j])), R.zero)));
}
