// Greatest common divisors of polynomials in several variables over a field.
//
// The method is the classical one: pick a variable, regard both polynomials
// as polynomials in that variable whose coefficients are polynomials in the
// rest, split off the contents (gcds of those coefficients, found by the same
// method with one variable fewer), and run Euclid's algorithm with
// pseudo-division on the primitive parts, taking the primitive part of each
// remainder so that coefficients stay small.

import { MPoly } from "./mpoly.js";
import { keyOf } from "./monomial.js";

const usedVars = (f, used) => {
  for (const { e } of f._m.values()) for (let i = 0; i < e.length; i += 1) if (e[i] > 0) used[i] = true;
};

// build Σ coeffs[k]·v^k
function fromCoefficients(R, i, coeffs) {
  const m = new Map();
  coeffs.forEach((c, k) => {
    for (const t of c._m.values()) {
      const e = [...t.e];
      e[i] += k;
      m.set(keyOf(e), { e, c: t.c });
    }
  });
  return new MPoly(R, m);
}

// pseudo-remainder of a by b as polynomials in variable i (b not free of i)
function pseudoRemainder(R, a, b, i) {
  const bc = b.coefficients(i);
  const db = bc.length - 1;
  const lb = bc[db];
  let r = a;
  for (;;) {
    const rc = r.coefficients(i);
    const dr = rc.length - 1;
    if (r.isZero() || dr < db) return r;
    // r ← lb·r − lc(r)·v^(dr−db)·b, which cancels the top coefficient
    const shift = new Array(R.n).fill(0);
    shift[i] = dr - db;
    const lead = rc[dr];
    const t = lead._m.size === 1 ? null : lead;
    let sub;
    if (t === null) {
      const [{ e, c }] = lead._m.values();
      const ee = e.map((x, k) => x + shift[k]);
      sub = b._mulTerm(ee, c);
    } else sub = b.mul(lead)._mulTerm(shift, R.K.one);
    r = r.mul(lb).sub(sub);
  }
}

function gcdList(R, list) {
  let g = R.zero;
  for (const c of list) {
    if (c.isZero()) continue;
    g = g.isZero() ? c : gcdPair(R, g, c);
    if (g.isConstant()) return R.one;
  }
  return g;
}

/** The content with respect to a variable: the gcd of the coefficients of its powers. */
function contentIn(R, f, i) {
  return gcdList(R, f.coefficients(i));
}

function gcdPair(R, f, g) {
  if (f.isZero()) return g.monic();
  if (g.isZero()) return f.monic();
  if (f.isConstant() || g.isConstant()) return R.one;
  // a monomial factor is quick to take out, and common in practice
  const used = new Array(R.n).fill(false);
  usedVars(f, used);
  const usedG = new Array(R.n).fill(false);
  usedVars(g, usedG);
  let i = -1;
  for (let k = 0; k < R.n; k += 1) {
    if (used[k] && usedG[k]) {
      i = k;
      break;
    }
  }
  if (i < 0) {
    // no shared variable: only the contents with respect to each other's variables can be common
    const k = used.findIndex(Boolean);
    return gcdPair(R, contentIn(R, f, k), g);
  }
  // a variable of f that g lacks: g can only share f's content in it
  for (let k = 0; k < R.n; k += 1) {
    if (used[k] && !usedG[k]) return gcdPair(R, contentIn(R, f, k), g);
    if (usedG[k] && !used[k]) return gcdPair(R, f, contentIn(R, g, k));
  }
  const cf = contentIn(R, f, i);
  const cg = contentIn(R, g, i);
  const c = gcdPair(R, cf, cg);
  let a = f.div(cf);
  let b = g.div(cg);
  if (a.degree(i) < b.degree(i)) [a, b] = [b, a];
  while (!b.isZero()) {
    const r = pseudoRemainder(R, a, b, i);
    a = b;
    b = r.isZero() ? r : r.div(contentIn(R, r, i));
  }
  // a is primitive in variable i; if it lost the variable, the primitive parts were coprime
  return (a.degree(i) > 0 ? c.mul(a) : c).monic();
}

/**
 * The greatest common divisor of two or more polynomials, normalised so its
 * leading coefficient is 1. gcd(0, 0) is 0.
 */
export function gcd(first, ...rest) {
  const list = Array.isArray(first) ? first : [first, ...rest];
  const base = list.find((f) => f instanceof MPoly);
  if (!base) throw new TypeError("gcd needs at least one polynomial");
  const R = base.ring;
  let g = R.zero;
  for (const f of list) g = gcdPair(R, g, R.from(f));
  return g;
}

/** The least common multiple, with leading coefficient 1. */
export function lcm(first, ...rest) {
  const list = Array.isArray(first) ? first : [first, ...rest];
  const base = list.find((f) => f instanceof MPoly);
  if (!base) throw new TypeError("lcm needs at least one polynomial");
  const R = base.ring;
  let l = R.one;
  for (const x of list) {
    const f = R.from(x);
    if (f.isZero()) return R.zero;
    l = l.mul(f).div(gcdPair(R, l, f)).monic();
  }
  return l;
}

/**
 * The content of f with respect to a variable (the gcd of the coefficients of
 * its powers, a polynomial in the other variables) and the primitive part
 * f / content.
 */
export function content(f, v) {
  const R = f.ring;
  return f.isZero() ? R.zero : contentIn(R, f, R.indexOf(v));
}

export function primitivePart(f, v) {
  return f.isZero() ? f : f.div(content(f, v));
}

/**
 * The square-free decomposition f = c·f₁·f₂²·f₃³…, with the fₖ square-free and
 * pairwise coprime: [{ factor, multiplicity }] and the constant c. Needs
 * characteristic 0 (derivatives vanish too easily otherwise).
 */
export function squareFree(f) {
  const R = f.ring;
  if (R.K.characteristic !== 0n && R.K.characteristic !== 0) throw new RangeError("squareFree needs a field of characteristic 0");
  if (f.isZero()) throw new RangeError("the zero polynomial has no factorisation");
  const unit = f._lead().c;
  let parts = [{ factor: f.monic(), multiplicity: 1 }];
  // split by one variable at a time with Yun's algorithm; factors free of that variable pass through
  for (let i = 0; i < R.n; i += 1) {
    const next = [];
    for (const { factor: p, multiplicity: m } of parts) {
      if (p.degree(i) <= 0) {
        next.push({ factor: p, multiplicity: m });
        continue;
      }
      const cont = contentIn(R, p, i);
      if (!cont.isConstant()) next.push({ factor: cont.monic(), multiplicity: m });
      const pp = p.div(cont);
      const d = pp.diff(i);
      let g = gcdPair(R, pp, d);
      let w = pp.div(g);
      let y = d.div(g);
      let k = 1;
      for (;;) {
        const z = y.sub(w.diff(i));
        if (z.isZero()) {
          if (!w.isConstant()) next.push({ factor: w.monic(), multiplicity: m * k });
          break;
        }
        const h = gcdPair(R, w, z);
        if (!h.isConstant()) next.push({ factor: h.monic(), multiplicity: m * k });
        w = w.div(h);
        y = z.div(h);
        k += 1;
      }
    }
    parts = next;
  }
  // merge equal multiplicities are kept separate only if the factors differ; combine same multiplicity
  const byMult = new Map();
  for (const { factor, multiplicity } of parts) {
    if (factor.isConstant()) continue;
    byMult.set(multiplicity, byMult.has(multiplicity) ? byMult.get(multiplicity).mul(factor) : factor);
  }
  const factors = [...byMult.entries()].sort((a, b) => a[0] - b[0]).map(([multiplicity, factor]) => ({ factor, multiplicity }));
  return { unit: R.K.toPublic(unit), factors };
}

/** The product of the distinct irreducible factors (characteristic 0). */
export function squareFreePart(f) {
  return squareFree(f).factors.reduce((p, { factor }) => p.mul(factor), f.ring.one);
}

export { fromCoefficients, pseudoRemainder };
