// Modules over a principal ideal domain, made concrete: a square matrix A
// over Q turns Qⁿ into a module over Q[x], with x acting as A. The structure
// theorem says Qⁿ ≅ Q[x]/(f₁) ⊕ … ⊕ Q[x]/(f_k) with f₁ | f₂ | … | f_k, and
// the fᵢ are the invariant factors of A: the Smith form of xI − A over Q[x].
// Each summand is cyclic, x acting on it by the companion matrix of fᵢ, so
// the blocks together are the rational canonical form. The elementary
// divisors (the prime-power pieces of the fᵢ) give the primary
// decomposition, and over an algebraically closed field they would be the
// Jordan blocks.

import { Rational } from "@ratmath/core";
import { smithForm, polynomialDomain } from "./smith.js";
import { factor } from "./factor.js";
import { fromRelations } from "./abelian.js";

const toQ = (x) => (x instanceof Rational ? x : new Rational(String(x)));

/** The companion matrix of a monic polynomial: ones below the diagonal, minus the coefficients in the last column. */
export function companionMatrix(f) {
  const d = f.degree();
  if (d < 1) return [];
  const c = f.coefficients(f.ring.vars[0]).map((p) => (p.isZero() ? new Rational(0n) : p.toConstant())); // constant first
  const M = Array.from({ length: d }, () => Array.from({ length: d }, () => new Rational(0n)));
  for (let i = 1; i < d; i += 1) M[i][i - 1] = new Rational(1n);
  for (let i = 0; i < d; i += 1) M[i][d - 1] = (c[i] ?? new Rational(0n)).negate();
  return M;
}

/**
 * The Q[x]-module of a square matrix A over Q (entries numbers, BigInts,
 * Core Rationals or text such as "1/2"). Returns
 * { invariantFactors, elementaryDivisors, minimalPolynomial,
 *   characteristicPolynomial, rationalCanonical, isCyclic, smith }.
 */
export function matrixModule(A, { variable = "x" } = {}) {
  const n = A.length;
  if (!n || A.some((r) => r.length !== n)) throw new RangeError("a square matrix is needed");
  const D = polynomialDomain(variable);
  const R = D.ring;
  const x = R.gen(0);
  const M = A.map((row, i) => row.map((a, j) => (i === j ? x.sub(R.from(toQ(a))) : R.from(toQ(a)).neg())));
  const S = smithForm(M, D);
  const inv = S.diagonal.filter((f) => f.degree() > 0);
  const elementary = [];
  for (const f of inv) {
    for (const { factor: g, multiplicity } of factor(f).factors) elementary.push({ factor: g.monic(), multiplicity, divisor: g.monic().pow(multiplicity) });
  }
  elementary.sort((a, b) => a.divisor.degree() - b.divisor.degree() || String(a.factor).localeCompare(String(b.factor)));
  const blocks = inv.map(companionMatrix);
  const C = Array.from({ length: n }, () => Array.from({ length: n }, () => new Rational(0n)));
  let at = 0;
  for (const B of blocks) {
    B.forEach((row, i) => row.forEach((v, j) => { C[at + i][at + j] = v; }));
    at += B.length;
  }
  const charpoly = S.diagonal.reduce((p, f) => p.mul(f), R.one);
  return Object.freeze({
    ring: R,
    invariantFactors: inv,
    elementaryDivisors: elementary,
    minimalPolynomial: inv.length ? inv[inv.length - 1] : R.one,
    characteristicPolynomial: charpoly,
    rationalCanonical: C,
    blocks,
    isCyclic: inv.length <= 1,
    smith: S,
    toString: () => (inv.length ? inv.map((f) => `Q[${variable}]/(${f})`).join(" ⊕ ") : "0"),
  });
}

/**
 * A finitely generated Z-module (abelian group) from generators and relations:
 * the quotient Zⁿ/L, L spanned by the rows. The same as fromRelations, under
 * the module name.
 */
export const integerModule = (relations, n) => fromRelations(relations, n);
