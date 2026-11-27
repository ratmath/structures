// Conway polynomials: the standard choice of an irreducible polynomial for
// each finite field GF(pⁿ), made so that the fields fit inside one another.
//
// The Conway polynomial C(p, n) is the first monic polynomial of degree n
// over F_p, in a fixed order, that is
//   - primitive: its root x generates the multiplicative group of GF(pⁿ), and
//   - compatible: for every m dividing n, x^((pⁿ − 1)/(p^m − 1)) is a root of
//     C(p, m), so the standard generator of the small field is a power of the
//     standard generator of the big one.
// The order compares xⁿ − a₍ₙ₋₁₎xⁿ⁻¹ + a₍ₙ₋₂₎xⁿ⁻² − … by the lists
// (a₍ₙ₋₁₎, …, a₀) with 0 < 1 < … < p − 1. Nothing here is looked up: the
// polynomials are found by search, which is practical while pⁿ is modest.

import { galoisField, isProbablePrime, primitiveRoot, factorInteger } from "@ratmath/algebraic";

const cache = new Map();

// polynomials over F_p as arrays of Numbers, constant first, with p < 2^20 so products stay exact
function mulmod(a, b, f, p) {
  const n = f.length - 1;
  const out = new Array(2 * n - 1).fill(0);
  for (let i = 0; i < a.length; i += 1) { if (!a[i]) continue; for (let j = 0; j < b.length; j += 1) out[i + j] = (out[i + j] + a[i] * b[j]) % p; }
  for (let i = out.length - 1; i >= n; i -= 1) {
    const c = out[i];
    if (!c) continue;
    for (let j = 0; j < n; j += 1) out[i - n + j] = (out[i - n + j] + (p - f[j]) * c) % p;
  }
  out.length = n;
  return out;
}
function powmod(e, f, p) {
  // x^e mod f
  const n = f.length - 1;
  let out = new Array(n).fill(0);
  out[0] = 1;
  if (n === 1) { // x ≡ −f[0]
    let r = 1n; let b = BigInt((p - f[0]) % p); const P = BigInt(p);
    for (let k = e; k > 0n; k >>= 1n) { if (k & 1n) r = (r * b) % P; b = (b * b) % P; }
    return [Number(r)];
  }
  let base = new Array(n).fill(0);
  base[1] = 1;
  for (let k = e; k > 0n; k >>= 1n) {
    if (k & 1n) out = mulmod(out, base, f, p);
    if (k > 1n) base = mulmod(base, base, f, p);
  }
  return out;
}
const isOne = (a) => a[0] === 1 && a.slice(1).every((c) => c === 0);

function primeFactors(n) {
  const out = [];
  for (const entry of factorInteger(n)) out.push(BigInt(Array.isArray(entry) ? entry[0] : entry.prime ?? entry.p ?? entry));
  return [...new Set(out.map(String))].map(BigInt);
}

/**
 * The Conway polynomial for GF(pⁿ): its coefficients, constant first, as
 * BigInts in 0 … p − 1 (the form galoisField takes as a modulus).
 */
export function conwayPolynomial(p, n, { limit = 2000000 } = {}) {
  p = Number(p);
  if (!Number.isInteger(p) || p < 2 || !isProbablePrime(BigInt(p))) throw new RangeError(`${p} is not prime`);
  if (!Number.isInteger(n) || n < 1) throw new RangeError("the degree is a positive whole number");
  if (p >= 1 << 20) throw new RangeError("Conway polynomials are computed here for primes below 2^20");
  const id = `${p}^${n}`;
  if (cache.has(id)) return [...cache.get(id)];
  const P = BigInt(p);
  const r = Number(primitiveRoot(P));
  let found;
  if (n === 1) found = [p - r, 1];
  else {
    if (p ** (n - 1) > limit) throw new RangeError(`finding the Conway polynomial for GF(${p}^${n}) would mean searching up to ${p}^${n - 1} candidates; the limit is ${limit}`);
    const N = P ** BigInt(n) - 1n;
    const cofactors = primeFactors(N).map((q) => N / q);
    const divisors = [];
    for (let q = 2; q <= n; q += 1) if (n % q === 0 && isProbablePrime(BigInt(q))) divisors.push(n / q);
    const smaller = divisors.map((m) => ({ m, C: conwayPolynomial(p, m, { limit }).map(Number), e: N / (P ** BigInt(m) - 1n) }));
    // a = (a₍ₙ₋₁₎, …, a₁) counts up in the Conway order; a₀ is forced: the norm of x is the generator for n = 1
    const a = new Array(n).fill(0);
    a[0] = r;
    const total = p ** (n - 1);
    for (let t = 0; t < total; t += 1) {
      let x = t;
      for (let i = 1; i < n; i += 1) { a[i] = x % p; x = Math.floor(x / p); }
      // coefficient of x^i is (−1)^(n−i)·aᵢ
      const f = a.map((c, i) => ((n - i) % 2 ? (p - c) % p : c));
      f.push(1);
      if (!isOne(powmod(N, f, p))) continue;
      if (cofactors.some((e) => isOne(powmod(e, f, p)))) continue;
      const compatible = smaller.every(({ C, e }) => {
        const y = powmod(e, f, p);
        let acc = new Array(n).fill(0);
        for (let i = C.length - 1; i >= 0; i -= 1) { acc = mulmod(acc, y, f, p); acc[0] = (acc[0] + C[i]) % p; }
        return acc.every((c) => c === 0);
      });
      if (compatible) { found = f; break; }
    }
    if (!found) throw new Error(`no Conway polynomial found for GF(${p}^${n})`);
  }
  const out = found.map(BigInt);
  cache.set(id, out);
  return [...out];
}

const SUP = "⁰¹²³⁴⁵⁶⁷⁸⁹";
const sup = (k) => String(k).split("").map((c) => SUP[Number(c)]).join("");

/** A polynomial over F_p (constant first) as text. */
export function formatModular(coeffs, v = "x") {
  const parts = [];
  for (let i = coeffs.length - 1; i >= 0; i -= 1) {
    const c = BigInt(coeffs[i]);
    if (c === 0n) continue;
    parts.push(i === 0 ? `${c}` : `${c === 1n ? "" : c}${v}${i > 1 ? sup(i) : ""}`);
  }
  return parts.join(" + ") || "0";
}

/**
 * GF(pⁿ) built on its Conway polynomial: a field object as from galoisField,
 * with `conway` (the polynomial), `text`, the standard generator `gen` (a
 * generator of the multiplicative group), and
 * `embed(m)` giving the map from the standard GF(p^m), m dividing n, into
 * this field: its standard generator goes to gen^((pⁿ − 1)/(p^m − 1)).
 */
export function conwayField(p, n, options) {
  const C = conwayPolynomial(p, n, options);
  const F = galoisField(p, n, n === 1 ? null : C);
  const P = BigInt(p);
  const gen = n === 1 ? (P - C[0]) % P : F.gen;
  const embed = (m) => {
    if (!Number.isInteger(m) || m < 1 || n % m) throw new RangeError(`GF(${p}^${m}) lies inside GF(${p}^${n}) only when ${m} divides ${n}`);
    const small = conwayField(p, m, options);
    const e = (P ** BigInt(n) - 1n) / (P ** BigInt(m) - 1n);
    const image = F.pow(gen, e);
    const map = (x) => {
      // x is a polynomial in the small generator (or a residue when m = 1)
      const coeffs = m === 1 ? [BigInt(x)] : x;
      let acc = F.zero;
      for (let i = coeffs.length - 1; i >= 0; i -= 1) acc = F.add(F.mul(acc, image), F.from(coeffs[i]));
      return acc;
    };
    return { from: small, exponent: e, image, map };
  };
  return { ...F, conway: C, text: formatModular(C), gen, generator: gen, embed };
}
