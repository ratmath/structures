// Factoring polynomials in several variables over the rationals.
//
// The method is Kronecker's. Substituting x₁ → t, x₂ → t^D, x₃ → t^(D²), …
// with D larger than every degree turns a polynomial in several variables
// into one in t without losing anything: each monomial lands on its own power
// of t, and the exponents can be read back as digits in base D. A factor of
// the original maps to a factor of the image, so: factor the image in one
// variable (which @ratmath/algebraic does), and try products of its
// irreducible factors, read back as polynomials in several variables, as
// divisors. A product of the fewest factors that divides is irreducible.
//
// The search is over subsets, so it is meant for the sizes one meets in
// teaching and in geometry, not for large inputs.

import { Rational } from "@ratmath/core";
import { factor as factorUnivariate } from "@ratmath/algebraic";
import { MPoly } from "./mpoly.js";
import { QQ } from "./field.js";
import { keyOf } from "./monomial.js";
import { squareFree } from "./gcd.js";

const gcdInt = (a, b) => {
  if (a < 0n) a = -a;
  if (b < 0n) b = -b;
  while (b) [a, b] = [b, a % b];
  return a;
};

/**
 * Write f = c·p with p having coprime integer coefficients and a positive
 * leading coefficient. Returns { content: c (a Rational), primitive: p }.
 */
export function integerPrimitive(f) {
  const R = f.ring;
  if (R.K !== QQ) throw new TypeError("integerPrimitive needs rational coefficients");
  if (f.isZero()) return { content: new Rational(0n, 1n), primitive: f };
  let num = 0n;
  let den = 1n;
  for (const { c } of f._m.values()) {
    num = gcdInt(num, c[0]);
    den = (den * c[1]) / gcdInt(den, c[1]);
  }
  let content = [num, den];
  if (f._lead().c[0] < 0n) content = [-num, den];
  return { content: new Rational(content[0], content[1]), primitive: f._scale(QQ.inv(content)) };
}

const MAX_IMAGE_DEGREE = 60000;
const MAX_PIECES = 20;

function factorSquareFree(p) {
  const R = p.ring;
  const used = p.variables();
  if (used.length === 0) return [];
  if (used.length === 1) {
    const i = R.indexOf(used[0]);
    const coeffs = p.coefficients(i).map((c) => (c.isZero() ? 0n : c._lead().c[0]));
    return factorUnivariate(coeffs).factors.map(({ factor: g }) => fromImage(R, g, [i], 0));
  }
  const idx = used.map((v) => R.indexOf(v));
  const D = Math.max(...idx.map((i) => p.degree(i))) + 1;
  const top = D ** idx.length;
  if (top > MAX_IMAGE_DEGREE) throw new RangeError(`${p} is too large to factor by Kronecker substitution (degrees up to ${D - 1} in ${idx.length} variables)`);
  const image = new Array(top).fill(0n);
  for (const { e, c } of p._m.values()) {
    let k = 0;
    idx.forEach((i, pos) => (k += e[i] * D ** pos));
    image[k] = c[0];
  }
  while (image.length && image[image.length - 1] === 0n) image.pop();
  const pieces = [];
  for (const { factor: g, multiplicity } of factorUnivariate(image).factors) for (let m = 0; m < multiplicity; m += 1) pieces.push(g);
  if (pieces.length > MAX_PIECES) throw new RangeError(`${p} has too many candidate factors (${pieces.length}) to search`);
  const found = [];
  let rest = p;
  let pool = pieces;
  let size = 1;
  while (2 * size <= pool.length) {
    let hit = null;
    for (const subset of subsets(pool.length, size)) {
      let prod = [1n];
      for (const k of subset) prod = mulInt(prod, pool[k]);
      const cand = fromImage(R, prod, idx, D);
      if (cand === null || cand.isConstant()) continue;
      const { quotient, remainder } = rest.divmod(cand);
      if (remainder.isZero()) {
        hit = { subset, cand, quotient };
        break;
      }
    }
    if (!hit) {
      size += 1;
      continue;
    }
    found.push(integerPrimitive(hit.cand).primitive);
    rest = integerPrimitive(hit.quotient).primitive;
    const gone = new Set(hit.subset);
    pool = pool.filter((_, k) => !gone.has(k));
  }
  if (!rest.isConstant()) found.push(rest);
  return found;
}

function mulInt(a, b) {
  const out = new Array(a.length + b.length - 1).fill(0n);
  for (let i = 0; i < a.length; i += 1) {
    if (a[i] === 0n) continue;
    for (let j = 0; j < b.length; j += 1) out[i + j] += a[i] * b[j];
  }
  return out;
}

// read a polynomial in t back through the substitution; D = 0 means one variable
function fromImage(R, coeffs, idx, D) {
  const m = new Map();
  for (let k = 0; k < coeffs.length; k += 1) {
    if (coeffs[k] === 0n) continue;
    const e = new Array(R.n).fill(0);
    if (D === 0) e[idx[0]] = k;
    else {
      let rest = k;
      for (let pos = 0; pos < idx.length; pos += 1) {
        e[idx[pos]] = rest % D;
        rest = Math.floor(rest / D);
      }
      if (rest !== 0) return null;
    }
    m.set(keyOf(e), { e, c: [coeffs[k], 1n] });
  }
  return new MPoly(R, m);
}

function* subsets(n, k, start = 0, acc = []) {
  if (acc.length === k) {
    yield acc;
    return;
  }
  for (let i = start; i <= n - (k - acc.length); i += 1) yield* subsets(n, k, i + 1, [...acc, i]);
}

/**
 * Factor a polynomial with rational coefficients into irreducibles over Q:
 *
 *   f = unit · Π factorᵢ ^ multiplicityᵢ
 *
 * Each factor has coprime integer coefficients and a positive leading
 * coefficient; `unit` is a Rational.
 */
export function factor(f) {
  const R = f.ring;
  if (R.K !== QQ) throw new TypeError(`factor works over the rationals; this polynomial is in ${R}`);
  if (f.isZero()) throw new RangeError("the zero polynomial has no factorisation");
  const { content } = integerPrimitive(f);
  const factors = [];
  for (const { factor: part, multiplicity } of squareFree(f).factors) {
    for (const g of factorSquareFree(integerPrimitive(part).primitive)) factors.push({ factor: g, multiplicity });
  }
  const cmp = R.cmp;
  factors.sort((a, b) => a.factor.degree() - b.factor.degree() || a.factor.length - b.factor.length || cmp(a.factor._lead().e, b.factor._lead().e) || (a.factor.toText() < b.factor.toText() ? -1 : 1));
  return { unit: content, factors };
}

/** Whether a non-constant polynomial over Q is irreducible. */
export function isIrreducible(f) {
  if (f.isConstant()) return false;
  const { factors } = factor(f);
  return factors.length === 1 && factors[0].multiplicity === 1;
}
