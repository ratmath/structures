// Determinants of matrices of polynomials, resultants and discriminants.
//
// The determinant is found by fraction-free (Bareiss) elimination: every
// division it performs is exact in the polynomial ring, so no fractions of
// polynomials ever appear and the entries stay as small as they can be.

import { MPoly } from "./mpoly.js";

/** The determinant of a square matrix whose entries are polynomials (or anything the ring reads). */
export function det(matrix, ringOf) {
  const R = ringOf ?? matrix.flat().find((x) => x instanceof MPoly)?.ring;
  if (!R) throw new TypeError("det needs a ring: pass one, or use polynomial entries");
  const n = matrix.length;
  const m = matrix.map((row) => {
    if (row.length !== n) throw new RangeError("the determinant needs a square matrix");
    return row.map((x) => R.from(x));
  });
  if (n === 0) return R.one;
  let sign = 1;
  let prev = R.one;
  for (let k = 0; k < n - 1; k += 1) {
    if (m[k][k].isZero()) {
      const swap = m.findIndex((row, i) => i > k && !row[k].isZero());
      if (swap < 0) return R.zero;
      [m[k], m[swap]] = [m[swap], m[k]];
      sign = -sign;
    }
    for (let i = k + 1; i < n; i += 1) {
      for (let j = k + 1; j < n; j += 1) {
        const t = m[i][j].mul(m[k][k]).sub(m[i][k].mul(m[k][j]));
        m[i][j] = k === 0 ? t : t.div(prev);
      }
      m[i][k] = R.zero;
    }
    prev = m[k][k];
  }
  const d = m[n - 1][n - 1];
  return sign < 0 ? d.neg() : d;
}

/**
 * The Sylvester matrix of f and g with respect to a variable: its
 * determinant is the resultant.
 */
export function sylvester(f, g, v) {
  const R = f.ring;
  g = R.from(g);
  const a = f.coefficients(v);
  const b = g.coefficients(v);
  if (!a.length || !b.length) throw new RangeError("the Sylvester matrix of the zero polynomial is not defined");
  const m = a.length - 1;
  const n = b.length - 1;
  const size = m + n;
  const rows = [];
  for (let i = 0; i < n; i += 1) rows.push(Array.from({ length: size }, (_, j) => (j >= i && j - i <= m ? a[m - (j - i)] : R.zero)));
  for (let i = 0; i < m; i += 1) rows.push(Array.from({ length: size }, (_, j) => (j >= i && j - i <= n ? b[n - (j - i)] : R.zero)));
  return rows;
}

/**
 * The resultant of f and g with respect to a variable: a polynomial in the
 * other variables that vanishes exactly where f and g have a common root in
 * that variable (or both leading coefficients vanish). It eliminates the
 * variable.
 */
export function resultant(f, g, v) {
  const R = f.ring;
  g = R.from(g);
  if (f.isZero() || g.isZero()) return R.zero;
  const i = R.indexOf(v);
  const m = f.degree(i);
  const n = g.degree(i);
  if (m === 0) return f.pow(n);
  if (n === 0) return g.pow(m);
  return det(sylvester(f, g, i), R);
}

/**
 * The discriminant with respect to a variable:
 * (−1)^(n(n−1)/2) · Res(f, ∂f/∂v) / lc(f). It vanishes exactly where f has a
 * repeated root in that variable.
 */
export function discriminant(f, v) {
  const R = f.ring;
  const i = R.indexOf(v);
  const n = f.degree(i);
  if (n < 1) throw new RangeError(`the discriminant needs degree at least 1 in ${R.vars[i]}`);
  const c = f.coefficients(i);
  const res = resultant(f, f.diff(i), i).div(c[n]);
  return ((n * (n - 1)) / 2) % 2 === 1 ? res.neg() : res;
}

/** The Jacobian matrix of a list of polynomials: rows are the gradients. */
export function jacobian(list) {
  return list.map((f) => f.gradient());
}

/** The Hessian matrix of second partial derivatives. */
export function hessian(f) {
  return f.ring.vars.map((u) => f.ring.vars.map((w) => f.diff(u).diff(w)));
}
