// Solving by radicals, and constructions with ruler and compass.
//
// A polynomial over Q can be solved by radicals exactly when its Galois
// group is solvable. When it is, this file writes the roots down: it walks a
// chain G = G₀ ⊳ G₁ ⊳ … ⊳ 1 with every step of prime index p, and at each
// step uses Lagrange's resolvents
//     rⱼ = Σₖ ωʲᵏ σᵏ(α),   ω a primitive p-th root of unity,
// for which rⱼᵖ lies one step down and α = (r₀ + r₁ + … + r₍ₚ₋₁₎)/p.
//
// Elements of the splitting field are never written symbolically. Each is
// held by its conjugates: the list of its certified complex values under
// every element of the Galois group, on which sums, products and the group
// act entrywise. Two facts turn those enclosures into exact answers:
//   - an element fixed by the whole group is rational, and a known multiple
//     of it is an integer, so a narrow enough enclosure names it;
//   - a non-zero element has a norm that is a non-zero integer over a known
//     denominator, so an element whose norm is certainly smaller is zero.
// ω is kept as a formal symbol (the algebra L ⊗ Q(ω)), so no assumption is
// made about which roots of unity the splitting field holds; ω itself is
// written in radicals by the same method applied to the cyclotomic field.
//
// The finished expression is then evaluated again from nothing but its
// rationals and certified n-th roots, and compared with the root.

import { Rational } from "@ratmath/core";
import { galoisAction, factor as factorZ } from "@ratmath/algebraic";
import { ring } from "./mpoly.js";

// ---- rationals as BigInt pairs

const gcdB = (a, b) => { a = a < 0n ? -a : a; b = b < 0n ? -b : b; while (b) [a, b] = [b, a % b]; return a; };
const normQ = (n, d) => { if (d < 0n) { n = -n; d = -d; } const g = gcdB(n, d) || 1n; return [n / g, d / g]; };
const Z = [0n, 1n];
const ONE = [1n, 1n];
const addQ = (a, b) => (a[0] === 0n ? b : b[0] === 0n ? a : normQ(a[0] * b[1] + b[0] * a[1], a[1] * b[1]));
const mulQ = (a, b) => (a[0] === 0n || b[0] === 0n ? Z : normQ(a[0] * b[0], a[1] * b[1]));
const absB = (a) => (a < 0n ? -a : a);
const lcmB = (a, b) => (a / gcdB(a, b)) * b;
const fdiv = (a, b) => { const q = a / b; return (a % b !== 0n && (a < 0n) !== (b < 0n)) ? q - 1n : q; };
const cdiv = (a, b) => -fdiv(-a, b);

class NeedPrecision extends Error {}

// ---- certified complex boxes [reLo, reHi, imLo, imHi] in fixed point at scale 2^P

function boxes(P) {
  const BP = BigInt(P);
  const one = 1n << BP;
  const up = (a) => -((-a) >> BP);
  const imul = (a0, a1, b0, b1) => {
    const p = [a0 * b0, a0 * b1, a1 * b0, a1 * b1];
    let lo = p[0];
    let hi = p[0];
    for (let i = 1; i < 4; i += 1) { if (p[i] < lo) lo = p[i]; if (p[i] > hi) hi = p[i]; }
    return [lo >> BP, up(hi)];
  };
  const add = (a, b) => [a[0] + b[0], a[1] + b[1], a[2] + b[2], a[3] + b[3]];
  const sub = (a, b) => [a[0] - b[1], a[1] - b[0], a[2] - b[3], a[3] - b[2]];
  const mul = (a, b) => {
    const rr = imul(a[0], a[1], b[0], b[1]);
    const ii = imul(a[2], a[3], b[2], b[3]);
    const ri = imul(a[0], a[1], b[2], b[3]);
    const ir = imul(a[2], a[3], b[0], b[1]);
    return [rr[0] - ii[1], rr[1] - ii[0], ri[0] + ir[0], ri[1] + ir[1]];
  };
  const scaleInt = (a, k) => (k >= 0n ? [a[0] * k, a[1] * k, a[2] * k, a[3] * k] : [a[1] * k, a[0] * k, a[3] * k, a[2] * k]);
  const scaleQ = (a, [n, d]) => { const b = scaleInt(a, n); return d === 1n ? b : [fdiv(b[0], d), cdiv(b[1], d), fdiv(b[2], d), cdiv(b[3], d)]; };
  const rational = ([n, d]) => { const lo = fdiv(n * one, d); const hi = cdiv(n * one, d); return [lo, hi, 0n, 0n]; };
  const interval = (lo, hi) => [fdiv(lo[0] * one, lo[1]), cdiv(hi[0] * one, hi[1])];
  const pow = (a, n) => { let out = rational(ONE); let b = a; for (let k = n; k > 0; k >>= 1) { if (k & 1) out = mul(out, b); if (k > 1) b = mul(b, b); } return out; };
  const maxAbs = (lo, hi) => (absB(lo) > absB(hi) ? absB(lo) : absB(hi));
  const minAbs = (lo, hi) => (lo <= 0n && hi >= 0n ? 0n : absB(lo) < absB(hi) ? absB(lo) : absB(hi));
  /** an upper bound and a lower bound for |z| over the box */
  const absUpper = (a) => maxAbs(a[0], a[1]) + maxAbs(a[2], a[3]);
  const absLower = (a) => { const r = minAbs(a[0], a[1]); const i = minAbs(a[2], a[3]); return r > i ? r : i; };
  const hasZero = (a) => a[0] <= 0n && a[1] >= 0n && a[2] <= 0n && a[3] >= 0n;
  const overlap = (a, b) => !(a[1] < b[0] || b[1] < a[0] || a[3] < b[2] || b[3] < a[2]);
  const toNumber = (x) => { const s = P > 60 ? BP - 60n : 0n; return Number(x >> s) / 2 ** Number(BP - s); };
  const mid = (a) => [toNumber((a[0] + a[1]) >> 1n), toNumber((a[2] + a[3]) >> 1n)];
  const fromNumber = (x) => { const k = BigInt(Math.round(x * 2 ** 50)); return P >= 50 ? k << (BP - 50n) : k >> (50n - BP); };
  const width = (a) => { const w = a[1] - a[0]; const h = a[3] - a[2]; return w > h ? w : h; };

  /** The n-th root of the box R nearest the complex number `near`, enclosed. */
  const rootNear = (R, n, near) => {
    let zr = fromNumber(near[0]);
    let zi = fromNumber(near[1]);
    const cr = (R[0] + R[1]) >> 1n;
    const ci = (R[2] + R[3]) >> 1n;
    const cmul = (a, b, c, d) => [(a * c - b * d) >> BP, (a * d + b * c) >> BP];
    const cpow = (a, b, k) => { let o = [one, 0n]; let q = [a, b]; for (let e = k; e > 0; e >>= 1) { if (e & 1) o = cmul(o[0], o[1], q[0], q[1]); if (e > 1) q = cmul(q[0], q[1], q[0], q[1]); } return o; };
    const steps = Math.ceil(Math.log2(P / 30)) + 3;
    for (let s = 0; s < steps; s += 1) {
      // z ← z − (zⁿ − c)/(n·zⁿ⁻¹)
      const [pr, pi] = cpow(zr, zi, n - 1);
      const [fr, fi] = cmul(pr, pi, zr, zi);
      const nr = fr - cr;
      const ni = fi - ci;
      const dr = pr * BigInt(n);
      const di = pi * BigInt(n);
      const den = (dr * dr + di * di) >> BP;
      if (den === 0n) throw new NeedPrecision();
      zr -= (nr * dr + ni * di) / den;
      zi -= (ni * dr - nr * di) / den;
    }
    const z = [zr, zr, zi, zi];
    const E = pow(z, n);
    const num = absUpper(sub(R, E));
    const low = absLower(E);
    if (low === 0n) throw new NeedPrecision();
    const h = (num * one) / low + 1n;
    if (h >= one >> 2n) throw new NeedPrecision();
    // the root is z(1 + δ) with |δ| ≤ 4h/n
    const rad = (absUpper(z) * 4n * h) / (one * BigInt(n)) + 2n;
    return [zr - rad, zr + rad, zi - rad, zi + rad];
  };
  return { P, one, add, sub, mul, scaleInt, scaleQ, rational, interval, pow, absUpper, absLower, hasZero, overlap, mid, width, rootNear, fromNumber };
}

// ---- elements of a Galois extension held by their conjugates

let nodeCount = 0;
/** A node of the expression graph: kind, operands, and `den` with den·(element) an algebraic integer. */
const node = (kind, den, a = null, b = null, extra = null) => ({ id: (nodeCount += 1), kind, den, a, b, extra, P: 0, val: null });

/**
 * A finite Galois extension of Q with its group listed: `size`, `mul[g][h]`
 * (the index of g∘h), `base(k, F)` the certified values of generator k under
 * every element, and a chain of subgroups.
 */
class Context {
  constructor({ size, mul, generators, chain, label }) {
    this.size = size;
    this.mul = mul;
    this.generators = generators; // [{ den, values(F) → boxes indexed by group element }]
    this.chain = chain; // [{ p, sigma, subgroup (set of indices of the smaller group) }], from the top down
    this.label = label;
    this.F = null;
    this.rationalCache = new Map();
    this.zeroCache = new Map();
    this.memo = new Map();
  }
  setPrecision(P) { this.F = boxes(P); }
  gen(k) { return node("gen", this.generators[k].den, k); }
  int(n) { return node("int", 1n, BigInt(n)); }
  add(x, y) { return node("add", lcmB(x.den, y.den), x, y); }
  sub(x, y) { return node("sub", lcmB(x.den, y.den), x, y); }
  mulN(x, y) { return node("mul", x.den * y.den, x, y); }
  /** the image of x under the group element h */
  act(h, x) { return h === 0 ? x : node("act", x.den, x, null, h); }

  /** The conjugates of x: its values under every group element, at the current precision. */
  values(x) {
    const F = this.F;
    if (x.P === F.P) return x.val;
    let out;
    const n = this.size;
    switch (x.kind) {
      case "gen": out = this.generators[x.a].values(F); break;
      case "int": { const b = F.rational([x.a, 1n]); out = Array(n).fill(b); break; }
      case "add": { const A = this.values(x.a); const B = this.values(x.b); out = A.map((v, i) => F.add(v, B[i])); break; }
      case "sub": { const A = this.values(x.a); const B = this.values(x.b); out = A.map((v, i) => F.sub(v, B[i])); break; }
      case "mul": { const A = this.values(x.a); const B = this.values(x.b); out = A.map((v, i) => F.mul(v, B[i])); break; }
      case "act": { const A = this.values(x.a); const h = x.extra; out = A.map((_, g) => A[this.mul[g][h]]); break; }
      default: throw new Error(`unknown node ${x.kind}`);
    }
    x.P = F.P;
    x.val = out;
    return out;
  }

  /** Whether x is exactly zero, decided through its norm. */
  isZero(x) {
    if (this.zeroCache.has(x.id)) return this.zeroCache.get(x.id);
    const F = this.F;
    const vals = this.values(x);
    let answer;
    if (vals.some((v) => !F.hasZero(v))) answer = false;
    else {
      // |norm of den·x| ≤ Π den·|value|; a non-zero integer is at least 1
      let bound = F.one;
      for (const v of vals) {
        bound = (bound * (F.absUpper(v) * x.den)) >> BigInt(F.P);
        if (bound === 0n) break;
      }
      if (bound >= F.one) throw new NeedPrecision();
      answer = true;
    }
    this.zeroCache.set(x.id, answer);
    return answer;
  }

  /** The rational value of an element fixed by the whole group. */
  rational(x) {
    if (this.rationalCache.has(x.id)) return this.rationalCache.get(x.id);
    const F = this.F;
    const v = this.values(x)[0];
    const s = F.scaleInt(v, x.den);
    if (s[1] - s[0] >= F.one >> 1n) throw new NeedPrecision();
    const m = (s[0] + (F.one >> 1n)) >> BigInt(F.P);
    let hit = null;
    for (const c of [m, m + 1n]) if (c * F.one >= s[0] && c * F.one <= s[1]) hit = c;
    if (hit === null || s[2] > 0n || s[3] < 0n) throw new Error("radicals: an element that must be rational is not; please report this polynomial");
    const q = normQ(hit, x.den);
    this.rationalCache.set(x.id, q);
    return q;
  }

  /** A double-precision value of x under the identity. */
  approx(x) { return this.F.mid(this.values(x)[0]); }
}

// ---- radical expressions: polynomials with rational coefficients in named radicals

const keyOf = (m) => m.map(([id, e]) => `${id}^${e}`).join("*");
class RPoly {
  constructor(atoms, terms = new Map()) { this.A = atoms; this.t = terms; }
  get isZero() { return this.t.size === 0; }
  /** the rational value when there are no radicals in it, else null */
  get constant() { if (this.t.size === 0) return Z; if (this.t.size === 1 && this.t.has("")) return this.t.get("").c; return null; }
  add(o) {
    const t = new Map(this.t);
    for (const [k, { c, m }] of o.t) {
      const s = addQ(t.get(k)?.c ?? Z, c);
      if (s[0] === 0n) t.delete(k); else t.set(k, { c: s, m });
    }
    return new RPoly(this.A, t);
  }
  scale(q) {
    if (q[0] === 0n) return new RPoly(this.A);
    const t = new Map();
    for (const [k, { c, m }] of this.t) t.set(k, { c: mulQ(c, q), m });
    return new RPoly(this.A, t);
  }
  neg() { return this.scale([-1n, 1n]); }
  mul(o) {
    let out = new RPoly(this.A);
    for (const { c: c1, m: m1 } of this.t.values()) for (const { c: c2, m: m2 } of o.t.values()) {
      const e = new Map(m1);
      for (const [id, k] of m2) e.set(id, (e.get(id) ?? 0) + k);
      out = out.add(this.A.monomial(mulQ(c1, c2), [...e.entries()].sort((x, y) => x[0] - y[0])));
    }
    return out;
  }
  pow(k) { let out = this.A.rat(ONE); for (let i = 0; i < k; i += 1) out = out.mul(this); return out; }
  key() { return [...this.t.entries()].map(([k, { c }]) => `${c[0]}/${c[1]}:${k}`).sort().join("+"); }
  atoms() { const s = new Set(); for (const { m } of this.t.values()) for (const [id] of m) s.add(id); return s; }
  /** q·P with P's coefficients whole, coprime, and its first term positive: [q, P] */
  content() {
    let num = 0n;
    let den = 1n;
    for (const { c } of this.t.values()) { num = gcdB(num, c[0]); den = lcmB(den, c[1]); }
    const first = [...this.t.keys()].sort()[0];
    const sign = this.t.get(first).c[0] < 0n ? -1n : 1n;
    const q = [sign * num, den];
    return [q, this.scale([q[1], q[0]])];
  }
}

const TWO_PI = 2 * Math.PI;
const SUP = "⁰¹²³⁴⁵⁶⁷⁸⁹";
const SUB = "₀₁₂₃₄₅₆₇₈₉";
const sup = (k) => String(k).split("").map((c) => SUP[Number(c)]).join("");
const sub = (k) => String(k).split("").map((c) => SUB[Number(c)]).join("");
const radicalSign = (n) => (n === 2 ? "√" : n === 3 ? "∛" : n === 4 ? "∜" : `${sup(n)}√`);
const fmtQ = ([n, d]) => (d === 1n ? `${absB(n)}` : `${absB(n)}/${d}`);

/** The radicals of one solution: each is an n-th root of an expression in the earlier ones, or a name for one. */
class Atoms {
  constructor() { this.reset(); }
  reset() { this.list = []; this.byKey = new Map(); this.unity = new Map(); this.rootCount = 0; }
  mark() { return this.list.length; }
  rewind(mark) {
    this.list.length = mark;
    for (const [key, ids] of this.byKey) this.byKey.set(key, ids.filter((id) => id < mark));
    for (const [p, v] of this.unity) if (v.mark >= mark) this.unity.delete(p);
    this.rootCount = this.list.filter((a) => a.name.startsWith("r")).length;
  }
  rat(q) { const p = new RPoly(this); if (q[0] !== 0n) p.t.set("", { c: q, m: [] }); return p; }
  atom(id, e = 1) { return this.monomial(ONE, [[id, e]]); }
  /** c·Π atomᵉ, with powers at or past an atom's index rewritten through what the atom is a root of */
  monomial(c, m) {
    for (let i = 0; i < m.length; i += 1) {
      const [id, e] = m[i];
      const a = this.list[id];
      if (a.kind === "unity" && e >= a.n) { const rest = m.map((x, j) => (j === i ? [id, e % a.n] : x)).filter(([, k]) => k > 0); return this.monomial(c, rest); }
      if (a.kind === "root" && e >= a.n) {
        const rest = m.map((x, j) => (j === i ? [id, e % a.n] : x)).filter(([, k]) => k > 0);
        return a.radicand.pow(Math.floor(e / a.n)).mul(this.monomial(c, rest));
      }
    }
    const p = new RPoly(this);
    const mm = m.filter(([, k]) => k > 0);
    if (c[0] !== 0n) p.t.set(keyOf(mm), { c, m: mm });
    return p;
  }
  /** The atom for the n-th root of `radicand` nearest `near`, or an earlier atom times a root of unity when that is the same number: [id, k] meaning atom·ζₙᵏ. */
  root(n, radicand, near) {
    const k = `${n}|${radicand.key()}`;
    for (const id of this.byKey.get(k) ?? []) {
      const a = this.list[id];
      const d = a.near[0] ** 2 + a.near[1] ** 2;
      const ratio = [(near[0] * a.near[0] + near[1] * a.near[1]) / d, (near[1] * a.near[0] - near[0] * a.near[1]) / d];
      const j = ((Math.round((Math.atan2(ratio[1], ratio[0]) * n) / TWO_PI) % n) + n) % n;
      return [id, j];
    }
    // a new radical stands for one definite root: the real one when the radicand is real and n is odd or
    // the radicand positive, otherwise the principal one (largest real part; of two, the one above the axis)
    const arg = Math.atan2(near[1], near[0]);
    const size = Math.hypot(near[0], near[1]);
    const argR = Math.atan2(Math.sin(arg * n), Math.cos(arg * n)); // the argument of the radicand, in (−π, π]
    let canon;
    if (Math.abs(Math.sin(argR)) < 1e-11 && (n % 2 === 1 || Math.cos(argR) > 0)) canon = Math.cos(argR) > 0 ? 0 : Math.PI;
    else canon = (Math.abs(Math.abs(argR) - Math.PI) < 1e-11 ? Math.PI : argR) / n;
    const j0 = ((Math.round(((arg - canon) * n) / TWO_PI) % n) + n) % n;
    near = [size * Math.cos(canon), size * Math.sin(canon)];
    if (Math.abs(near[1]) < 1e-13 * size) near[1] = 0;
    if (Math.abs(near[0]) < 1e-13 * size) near[0] = 0;
    const id = this.list.length;
    const q = radicand.constant;
    let name;
    if (q && q[1] === 1n && q[0] > 0n) name = `${radicalSign(n)}${q[0]}`;
    else if (q && q[0] === -1n && q[1] === 1n && n === 2) name = "i";
    else { this.rootCount += 1; name = `r${sub(this.rootCount)}`; }
    this.list.push({ id, kind: "root", n, radicand, near, name, simple: !name.startsWith("r") });
    this.byKey.set(k, [...(this.byKey.get(k) ?? []), id]);
    return [id, j0];
  }
  /** A named primitive p-th root of unity (its p-th power is 1), standing for `expr`. */
  unityAtom(p, expr) {
    const id = this.list.length;
    this.list.push({ id, kind: "unity", n: p, expr, near: [Math.cos(TWO_PI / p), Math.sin(TWO_PI / p)], name: `ζ${sub(p)}`, simple: false });
    return id;
  }
}

function formatRPoly(p, atoms, { tex = false, inline = false, budget = { left: Infinity } } = {}) {
  if (p.isZero) return "0";
  let L = 1n;
  for (const { c } of p.t.values()) L = lcmB(L, c[1]);
  const name = (id) => {
    const a = atoms.list[id];
    if (inline && !a.simple) {
      if (a.kind === "unity") return `(${formatRPoly(a.expr, atoms, { tex, inline, budget })})`;
      const inner = formatRPoly(a.radicand, atoms, { tex, inline, budget });
      return tex ? (a.n === 2 ? `\\sqrt{${inner}}` : `\\sqrt[${a.n}]{${inner}}`) : `${radicalSign(a.n)}(${inner})`;
    }
    if (!tex) return a.name;
    if (a.kind === "unity") return a.name.replace(/ζ([₀-₉]+)/, (_, d) => `\\zeta_{${[...d].map((ch) => SUB.indexOf(ch)).join("")}}`);
    if (a.simple) { if (a.name === "i") return "i"; const q = a.radicand.constant; return a.n === 2 ? `\\sqrt{${q[0]}}` : `\\sqrt[${a.n}]{${q[0]}}`; }
    return a.name.replace(/r([₀-₉]+)/, (_, d) => `r_{${[...d].map((ch) => SUB.indexOf(ch)).join("")}}`);
  };
  const terms = [...p.t.values()].sort((x, y) => x.m.length - y.m.length || keyOf(x.m).localeCompare(keyOf(y.m)));
  const parts = [];
  for (const { c, m } of terms) {
    const k = (c[0] * L) / c[1];
    const neg = k < 0n;
    // i first; the square roots of whole numbers in one monomial are written as one square root
    const isSqrt = ([id, e]) => { const a = atoms.list[id]; return e === 1 && a.kind === "root" && a.n === 2 && a.simple && a.name !== "i"; };
    const hasI = m.some(([id]) => atoms.list[id].name === "i");
    const merged = m.filter(isSqrt).reduce((acc, [id]) => acc * atoms.list[id].radicand.constant[0], 1n);
    const rest = m.filter((x) => !isSqrt(x) && atoms.list[x[0]].name !== "i").map(([id, e]) => `${name(id)}${e > 1 ? (tex ? `^{${e}}` : sup(e)) : ""}`);
    const pieces = [...(merged > 1n ? [tex ? `\\sqrt{${merged}}` : `√${merged}`] : []), ...rest];
    const mono = (hasI ? "i" : "") + (hasI && pieces.length && !tex ? "·" : hasI && pieces.length ? " " : "") + pieces.join(tex ? " " : "·");
    const coef = absB(k) === 1n && mono ? "" : `${absB(k)}`;
    const body = `${coef}${coef && mono && !tex ? "·" : coef && mono ? " " : ""}${mono}` || "1";
    parts.push((parts.length ? (neg ? (tex ? " - " : " − ") : " + ") : neg ? (tex ? "-" : "−") : "") + body);
    budget.left -= body.length + 3;
    if (budget.left < 0) throw new RangeError("too long");
  }
  const body = parts.join("");
  if (L === 1n) return body;
  if (tex) return `\\frac{${body}}{${L}}`;
  return parts.length === 1 ? `${body}/${L}` : `(${body})/${L}`;
}

// ---- the Lagrange resolvent recursion

/** n = sᵏ·t with t free of k-th powers, by trial division (enough to tidy small radicands). */
function extractPower(n, k) {
  let s = 1n;
  let t = n;
  for (let p = 2n; p < 2000n && p ** BigInt(k) <= t; p += p === 2n ? 1n : 2n) {
    const pk = p ** BigInt(k);
    while (t % pk === 0n) { t /= pk; s *= p; }
  }
  return [s, t];
}

class Solver {
  constructor(atoms) {
    this.atoms = atoms;
    this.unityContexts = new Map();
  }

  /** ζ_p for a prime p, as an expression in radicals (found once): written out for 3 and 5, named from 7 on. */
  unity(p) {
    const A = this.atoms;
    if (p === 2) return A.rat([-1n, 1n]);
    if (A.unity.has(p)) return A.unity.get(p).poly;
    const mark = A.mark();
    const ctx = cyclotomicContext(p);
    ctx.setPrecision(this.P);
    const expr = this.express(ctx, ctx.gen(0), ctx.chain.length);
    const poly = p <= 5 ? expr : A.atom(A.unityAtom(p, expr));
    A.unity.set(p, { poly, mark });
    return poly;
  }

  /** ζₙᵏ for prime n. */
  unityPower(n, k) {
    if (n === 2) return this.atoms.rat(k % 2 ? [-1n, 1n] : ONE);
    return this.unity(n).pow(((k % n) + n) % n);
  }

  /** √m for a whole number m ≠ 0, as (a whole number)·[i]·Π √p over primes: [coefficient, product of atoms] */
  sqrtInteger(m) {
    const A = this.atoms;
    let out = A.rat(ONE);
    if (m < 0n) out = out.mul(A.atom(A.root(2, A.rat([-1n, 1n]), [0, 1])[0]));
    let t = absB(m);
    let s = 1n;
    const take = (p) => { let e = 0; while (t % p === 0n) { t /= p; e += 1; } s *= p ** BigInt(e >> 1); if (e & 1) out = out.mul(A.atom(A.root(2, A.rat([p, 1n]), [Math.sqrt(Number(p)), 0])[0])); };
    for (let p = 2n; p * p <= t && p < 100000n; p += p === 2n ? 1n : 2n) if (t % p === 0n) take(p);
    if (t > 1n) take(t);
    return out.scale([s, 1n]);
  }

  /** The n-th root of an expression nearest `near` (n prime), tidied: whole powers pulled out, equal radicals shared. */
  radical(n, radicand, near) {
    const A = this.atoms;
    if (radicand.isZero) return new RPoly(A);
    const [q, P] = radicand.content();
    // q = m/d = (m·dⁿ⁻¹)/dⁿ
    const m = q[0] * q[1] ** BigInt(n - 1);
    const constant = P.constant; // then P = 1
    if (n === 2 && constant) {
      // ±√m/d with the sign that lands nearest
      const root = this.sqrtInteger(m).scale([1n, q[1]]);
      const positive = m > 0n ? near[0] >= 0 : near[1] >= 0;
      return positive ? root : root.neg();
    }
    const [s, t] = extractPower(absB(m), n);
    // ⁿ√(q·P) = (s/d)·ⁿ√(±t·P); for odd n the sign comes out as well
    const outside = n % 2 === 1 && m < 0n ? [-s, q[1]] : [s, q[1]];
    const inside = P.scale([n % 2 === 0 && m < 0n ? -t : t, 1n]);
    const f = Number(outside[0]) / Number(outside[1]);
    const [id, k] = A.root(n, inside, [near[0] / f, near[1] / f]);
    let out = A.atom(id).scale(outside);
    if (k !== 0) out = out.mul(this.unityPower(n, k));
    return out;
  }

  /**
   * For j = 1 … p−1, whether the number rⱼ = Σₖ ζᵏʲ·sₖ is non-zero. A non-zero
   * rⱼ is an algebraic integer over a known denominator of degree at most
   * |G|(p−1), with every conjugate among the values Σₖ ζᵏⁱ·g(sₖ); so it cannot
   * be smaller than the reciprocal of the product of the others.
   */
  liveResolvents(ctx, s, p, below) {
    this.rationalPowers = [];
    const F = ctx.F;
    const one = F.rational(ONE);
    const zeta = Array.from({ length: p }, (_, k) => (k === 0 ? one : F.rootNear(one, p, [Math.cos((TWO_PI * k) / p), Math.sin((TWO_PI * k) / p)])));
    const vals = s.map((v) => ctx.values(v));
    const den = s[0].den;
    const at = (g, j) => { let acc = F.rational(Z); for (let k = 0; k < p; k += 1) acc = F.add(acc, F.mul(zeta[(j * k) % p], vals[k][g])); return F.scaleInt(acc, den); };
    let M = F.one;
    for (let g = 0; g < ctx.size; g += 1) for (let j = 1; j < p; j += 1) { const u = F.absUpper(at(g, j)); if (u > M) M = u; }
    const degree = ctx.size * (p - 1);
    const out = [false];
    for (let j = 1; j < p; j += 1) {
      const y = at(0, j);
      if (!F.hasZero(y)) {
        out.push(true);
        // is rⱼᵖ a rational number? then (den·rⱼ)ᵖ is a whole number m, and the difference from m is tested as above
        const Y = F.pow(y, p);
        const m = (((Y[0] + Y[1]) >> 1n) + (F.one >> 1n)) >> BigInt(F.P);
        const diff = F.sub(Y, F.rational([m, 1n]));
        if (F.hasZero(diff)) {
          let big = F.one;
          for (let i = 0; i < p; i += 1) big = (big * M) >> BigInt(F.P);
          big += absB(m) * F.one;
          let bound = F.absUpper(diff);
          // rⱼᵖ lies in the field one step down with ζ adjoined, of degree at most `below`·(p − 1)
          for (let i = 1; i < below * (p - 1) && bound > 0n; i += 1) bound = (bound * big) >> BigInt(F.P);
          if (bound >= F.one) throw new NeedPrecision();
          this.rationalPowers[j] = normQ(m, den ** BigInt(p));
        }
        continue;
      }
      let bound = F.absUpper(y);
      for (let i = 1; i < degree && bound > 0n; i += 1) bound = (bound * M) >> BigInt(F.P);
      if (bound >= F.one) throw new NeedPrecision();
      out.push(false);
    }
    return out;
  }

  /** x, an element of the fixed field of the level-th subgroup of the chain, in radicals. */
  express(ctx, x, level) {
    if (level === 0) return this.atoms.rat(ctx.rational(x));
    const memoKey = `${x.id}@${level}`;
    if (ctx.memo.has(memoKey)) return ctx.memo.get(memoKey);
    const { p, sigma } = ctx.chain[level - 1];
    // σᵏ(x)
    const s = [x];
    let h = sigma;
    for (let k = 1; k < p; k += 1) { s.push(ctx.act(h, x)); h = ctx.mul[sigma][h]; }
    let out;
    if (ctx.isZero(ctx.sub(s[1], x))) out = this.express(ctx, x, level - 1); // already one step down
    else {
      const r0 = s.reduce((a, b) => ctx.add(a, b));
      out = this.express(ctx, r0, level - 1);
      if (p === 2) {
        const r1 = ctx.sub(s[0], s[1]);
        const R = this.express(ctx, ctx.mulN(r1, r1), level - 1);
        out = out.add(this.radical(2, R, ctx.approx(r1)));
      } else {
        // r₁ = Σ ωᵏ σᵏ(x) in L[ω]/(1 + ω + … + ω^(p−1)): coefficients of 1, ω, …, ω^(p−2)
        const reduceTop = (c) => c.slice(0, p - 1).map((v) => (c[p - 1] ? (v ? ctx.sub(v, c[p - 1]) : ctx.sub(ctx.int(0), c[p - 1])) : v ?? ctx.int(0)));
        const mulW = (A, B) => {
          const c = Array(p).fill(null);
          for (let i = 0; i < p - 1; i += 1) for (let j = 0; j < p - 1; j += 1) {
            const t = ctx.mulN(A[i], B[j]);
            const k = (i + j) % p;
            c[k] = c[k] ? ctx.add(c[k], t) : t;
          }
          return reduceTop(c);
        };
        const r1 = reduceTop(s);
        let R1 = null;
        let base = r1;
        for (let e = p; e > 0; e >>= 1) { if (e & 1) R1 = R1 ? mulW(R1, base) : base; if (e > 1) base = mulW(base, base); }
        // which resolvents are not zero as numbers (ω may satisfy relations with the field that the symbol does not)
        const live = this.liveResolvents(ctx, s, p, ctx.chain.slice(0, level - 1).reduce((a, c) => a * c.p, 1));
        const rationalPowers = this.rationalPowers;
        if (live.some(Boolean)) {
          // the coefficients of r₁ᵖ lie one step down
          const coeffs = R1.map((c) => (ctx.isZero(c) ? new RPoly(this.atoms) : this.express(ctx, c, level - 1)));
          const w = this.unity(p);
          const wPow = [this.atoms.rat(ONE)];
          for (let k = 1; k < p; k += 1) wPow.push(wPow[k - 1].mul(w));
          const sValues = s.map((v) => ctx.approx(v));
          for (let j = 1; j < p; j += 1) {
            if (!live[j]) continue;
            const known = rationalPowers[j];
            // rⱼ is r₁ with ω replaced by ωʲ, and so is its p-th power
            let Rj = new RPoly(this.atoms);
            if (known) Rj = this.atoms.rat(known);
            else coeffs.forEach((c, l) => { if (!c.isZero) Rj = Rj.add(c.mul(wPow[(j * l) % p])); });
            let re = 0;
            let im = 0;
            sValues.forEach(([a, b], k) => { const t = (TWO_PI * ((j * k) % p)) / p; re += a * Math.cos(t) - b * Math.sin(t); im += a * Math.sin(t) + b * Math.cos(t); });
            out = out.add(this.radical(p, Rj, [re, im]));
          }
        }
      }
      out = out.scale([1n, BigInt(p)]);
    }
    ctx.memo.set(memoKey, out);
    return out;
  }
}

// ---- the two kinds of Galois extension used

const primeFactorsSmall = (n) => { const out = []; for (let p = 2; p * p <= n; p += 1) while (n % p === 0) { out.push(p); n /= p; } if (n > 1) out.push(n); return out; };

/** Q(ζ_p) for an odd prime p, with group (Z/p)* and the chain of its subgroups, the largest prime steps first. */
function cyclotomicContext(p) {
  const size = p - 1;
  // elements a = 1 … p−1 at index a−1 (the identity first)
  const mul = Array.from({ length: size }, (_, i) => Array.from({ length: size }, (_, j) => (((i + 1) * (j + 1)) % p) - 1));
  let g = 2;
  const orderOf = (a) => { let k = 1; let x = a; while (x !== 1) { x = (x * a) % p; k += 1; } return k; };
  while (orderOf(g) !== size) g += 1;
  const primes = primeFactorsSmall(size).sort((a, b) => b - a);
  const chain = [];
  let index = 1; // G_level = ⟨g^index⟩
  for (const q of primes) {
    let sigma = 1;
    for (let k = 0; k < index; k += 1) sigma = (sigma * g) % p;
    chain.push({ p: q, sigma: sigma - 1 });
    index *= q;
  }
  const generators = [{
    den: 1n,
    values: (F) => {
      const one = F.rational(ONE);
      return Array.from({ length: size }, (_, i) => {
        const t = (TWO_PI * (i + 1)) / p;
        return F.rootNear(one, p, [Math.cos(t), Math.sin(t)]);
      });
    },
  }];
  return new Context({ size, mul, generators, chain, label: `Q(ζ${sub(p)})` });
}

/** A chain G ⊳ … ⊳ 1 with steps of prime index for a solvable group given by its multiplication table; null when not solvable. */
function primeChain(size, mul) {
  const inv = Array(size);
  for (let g = 0; g < size; g += 1) for (let h = 0; h < size; h += 1) if (mul[g][h] === 0) inv[g] = h;
  const closure = (gens, start = [0]) => {
    const seen = new Set(start);
    const queue = [...seen];
    const all = [...gens];
    for (let i = 0; i < queue.length; i += 1) for (const g of all) { const y = mul[g][queue[i]]; if (!seen.has(y)) { seen.add(y); queue.push(y); } }
    return seen;
  };
  const derived = (H) => {
    const els = [...H];
    const comms = new Set();
    for (const a of els) for (const b of els) comms.add(mul[mul[a][b]][mul[inv[a]][inv[b]]]);
    return closure([...comms]);
  };
  const steps = []; // from the top down: { p, sigma, group (the larger one) }
  let H = new Set(Array.from({ length: size }, (_, i) => i));
  while (H.size > 1) {
    const D = derived(H);
    if (D.size === H.size) return null;
    // climb from D to H by prime steps; each subgroup between is normal in H because H/D is abelian
    const up = [];
    let cur = D;
    while (cur.size < H.size) {
      const g = [...H].find((x) => !cur.has(x));
      // the order of g modulo cur
      let m = 1;
      let x = g;
      while (!cur.has(x)) { x = mul[g][x]; m += 1; }
      const q = primeFactorsSmall(m)[0];
      let sigma = 0;
      for (let k = 0; k < m / q; k += 1) sigma = mul[g][sigma];
      const next = closure([...cur, sigma], [...cur]);
      up.push({ p: q, sigma });
      cur = next;
    }
    steps.push(...up.reverse());
    H = D;
  }
  return steps;
}

function polynomialContext(action) {
  const els = action.group.elements;
  const n = action.degree;
  const keyP = (g) => g.join(",");
  const idKey = keyP(Array.from({ length: n }, (_, i) => i));
  const order = [els.findIndex((g) => keyP(g) === idKey), ...els.keys()].filter((i, j, arr) => j === 0 || i !== arr[0]);
  const perms = order.map((i) => els[i]);
  const index = new Map(perms.map((g, i) => [keyP(g), i]));
  const size = perms.length;
  const mul = perms.map((g) => perms.map((h) => index.get(keyP(h.map((x) => g[x])))));
  const chain = primeChain(size, mul);
  const lead = absB(action.polynomial[n]);
  const f = action.polynomial;
  // real roots are refined by @ratmath/algebraic; the others by Newton's method here, then enclosed:
  // a polynomial of degree n has a root within n·|f(z)/f′(z)| of any point z
  const rootBoxes = (F) => action.roots.map((x) => {
    if (x.rational) return F.rational(x.rational);
    if (!x.im) { x.refine(F.P + 8); const re = F.interval(x.lo, x.hi); return [re[0], re[1], 0n, 0n]; }
    x.refine(64);
    const coarse = [...F.interval(x.re[0], x.re[1]), ...F.interval(x.im[0], x.im[1])];
    const BP = BigInt(F.P);
    let zr = (coarse[0] + coarse[1]) >> 1n;
    let zi = (coarse[2] + coarse[3]) >> 1n;
    const steps = Math.ceil(Math.log2(F.P / 40)) + 3;
    for (let t = 0; t < steps; t += 1) {
      let [fr, fi, dr, di] = [0n, 0n, 0n, 0n];
      for (let i = n; i >= 0; i -= 1) {
        [dr, di] = [((dr * zr - di * zi) >> BP) + fr, ((dr * zi + di * zr) >> BP) + fi];
        [fr, fi] = [((fr * zr - fi * zi) >> BP) + f[i] * F.one, (fr * zi + fi * zr) >> BP];
      }
      const den = (dr * dr + di * di) >> BP;
      if (den === 0n) throw new NeedPrecision();
      zr -= (fr * dr + fi * di) / den;
      zi -= (fi * dr - fr * di) / den;
    }
    const z = [zr, zr, zi, zi];
    let fb = F.rational(Z);
    let db = F.rational(Z);
    for (let i = n; i >= 0; i -= 1) { db = F.add(F.mul(db, z), fb); fb = F.add(F.mul(fb, z), F.rational([f[i], 1n])); }
    const low = F.absLower(db);
    if (low === 0n) throw new NeedPrecision();
    const rad = (BigInt(n) * F.absUpper(fb) * F.one) / low + 2n;
    const box = [zr - rad, zr + rad, zi - rad, zi + rad];
    if (box[0] < coarse[0] || box[1] > coarse[1] || box[2] < coarse[2] || box[3] > coarse[3]) throw new NeedPrecision();
    return box;
  });
  let cache = { P: 0, boxes: null };
  const generators = Array.from({ length: n }, (_, k) => ({
    den: lead,
    values: (F) => {
      if (cache.P !== F.P) cache = { P: F.P, boxes: rootBoxes(F) };
      return perms.map((g) => cache.boxes[g[k]]);
    },
  }));
  return { ctx: chain ? new Context({ size, mul, generators, chain, label: "splitting field" }) : null, perms, rootBoxes };
}

// ---- results

const decimal = (x, digits) => {
  // x is a fixed-point BigInt at scale 2^P handed in as [value, P]
  const [v, P] = x;
  const neg = v < 0n;
  const a = neg ? -v : v;
  const scaled = (a * 10n ** BigInt(digits)) >> BigInt(P);
  const s = scaled.toString().padStart(digits + 1, "0");
  return `${neg ? "−" : ""}${s.slice(0, -digits)}.${s.slice(-digits)}`;
};

/** One number written in radicals. */
export class RadicalExpression {
  #atoms;
  #expr;
  #target;
  constructor(atoms, expr, target) {
    this.#atoms = atoms;
    this.#expr = expr;
    this.#target = target; // (F) → certified box of the number meant, or null
    Object.freeze(this);
  }

  #needed() {
    const need = new Set();
    const visit = (p) => { for (const id of p.atoms()) if (!need.has(id)) { need.add(id); const a = this.#atoms.list[id]; visit(a.kind === "unity" ? a.expr : a.radicand); } };
    visit(this.#expr);
    return [...need].sort((a, b) => a - b).map((id) => this.#atoms.list[id]);
  }

  /** The rational value when no radicals are needed, else null. */
  get rational() { const q = this.#expr.constant; return q ? new Rational(q[0], q[1]) : null; }

  /**
   * The radicals the expression is built from, innermost first:
   * [{ name, kind ("root" or "name"), index, of, value }]. A root is the
   * `index`-th root of `of` lying nearest `value`, which settles the branch.
   */
  get definitions() {
    return this.#needed().filter((a) => !a.simple).map((a) => ({
      name: a.name,
      kind: a.kind === "root" ? "root" : "name",
      index: a.kind === "root" ? a.n : null,
      of: formatRPoly(a.kind === "root" ? a.radicand : a.expr, this.#atoms),
      value: { re: a.near[0], im: a.near[1] },
      text: a.kind === "root" ? `${a.name} = ${radicalSign(a.n)}(${formatRPoly(a.radicand, this.#atoms)})` : `${a.name} = ${formatRPoly(a.expr, this.#atoms)}`,
    }));
  }

  /** The expression in terms of the named radicals of `definitions`. */
  get text() { return formatRPoly(this.#expr, this.#atoms); }

  /** The whole expression written out with nested radicals, or null when it would pass `maxLength` characters. */
  inline({ maxLength = 600, tex = false } = {}) {
    try { return formatRPoly(this.#expr, this.#atoms, { inline: true, tex, budget: { left: maxLength } }); } catch (e) { if (e instanceof RangeError) return null; throw e; }
  }
  toTeX(options = {}) { return this.inline({ ...options, tex: true }); }

  /** The largest index of a root used, and whether only square roots are. */
  get indices() { return [...new Set(this.#needed().filter((a) => a.kind === "root").map((a) => a.n))].sort((a, b) => a - b); }
  get usesOnlySquareRoots() { return this.indices.every((n) => n === 2); }
  /** How many radicals the expression uses. */
  get size() { return this.#needed().filter((a) => a.kind === "root").length; }

  #evaluate(F) {
    const val = new Map();
    const evalPoly = (p) => {
      let acc = F.rational(Z);
      for (const { c, m } of p.t.values()) {
        let t = F.rational(c);
        for (const [id, e] of m) t = F.mul(t, F.pow(val.get(id), e));
        acc = F.add(acc, t);
      }
      return acc;
    };
    for (const a of this.#needed()) val.set(a.id, a.kind === "unity" ? evalPoly(a.expr) : F.rootNear(evalPoly(a.radicand), a.n, a.near));
    return evalPoly(this.#expr);
  }

  /**
   * The value, computed from the rationals and certified roots alone:
   * { re, im } as decimal strings correct to the digits shown, give or take
   * one in the last place.
   */
  evaluate(digits = 30) {
    for (let P = Math.max(128, Math.ceil(digits * 3.33) + 64); P <= 1 << 15; P *= 2) {
      try {
        const F = boxes(P);
        const b = this.#evaluate(F);
        if (F.width(b) > F.one >> BigInt(Math.ceil(digits * 3.33) + 4)) continue;
        return { re: decimal([(b[0] + b[1]) >> 1n, P], digits), im: decimal([(b[2] + b[3]) >> 1n, P], digits), real: b[2] <= 0n && b[3] >= 0n };
      } catch (e) { if (!(e instanceof NeedPrecision)) throw e; }
    }
    throw new Error("radicals: could not evaluate the expression to that many digits");
  }

  /** A floating-point value: a number when the imaginary part is negligible, else { re, im }. */
  approx() {
    const v = this.evaluate(20);
    const re = Number(v.re.replace("−", "-"));
    const im = Number(v.im.replace("−", "-"));
    return Math.abs(im) < 1e-15 ? re : { re, im };
  }

  /**
   * Evaluates the expression afresh and compares it with the number it
   * stands for: { agrees, digits }. `digits` is how many decimal places the
   * two certified enclosures were compared to.
   */
  check(digits = 40) {
    if (!this.#target) return { agrees: null, digits: 0 };
    for (let P = Math.ceil(digits * 3.33) + 64; P <= 1 << 15; P *= 2) {
      try {
        const F = boxes(P);
        const mine = this.#evaluate(F);
        const { box, others } = this.#target(F);
        if (F.width(mine) > F.one >> BigInt(Math.ceil(digits * 3.33))) continue;
        // it must meet the root meant and miss every other root
        if (others.some((o) => F.overlap(mine, o))) continue;
        return { agrees: F.overlap(mine, box), digits };
      } catch (e) { if (!(e instanceof NeedPrecision)) throw e; }
    }
    return { agrees: null, digits: 0 };
  }

  toString() {
    const inline = this.inline({ maxLength: 160 });
    if (inline) return inline;
    return [...this.definitions.map((d) => d.text), `value = ${this.text}`].join("; ");
  }
}

function withPrecision(run) {
  for (let P = 192; P <= 1 << 15; P *= 2) {
    try { return run(P); } catch (e) { if (!(e instanceof NeedPrecision)) throw e; }
  }
  throw new Error("radicals: the precision limit was reached");
}

function readPolynomial(input) {
  if (typeof input === "string") {
    const letters = [...new Set(input.match(/[A-Za-z_Ͱ-Ͽ]\w*/g) ?? [])];
    if (letters.length !== 1) throw new RangeError("give a polynomial in one variable, such as \"x^3 - 2\"");
    const R = ring([letters[0]]);
    return R.poly(input).coefficients(letters[0]).map((c) => c.toConstant());
  }
  if (Array.isArray(input)) return input;
  if (input && typeof input.toUnivariate === "function") return input.toUnivariate();
  throw new TypeError("a polynomial is text in one variable or a list of coefficients, constant first");
}

function readQ(c) {
  if (Array.isArray(c)) return normQ(BigInt(c[0]), BigInt(c[1]));
  if (typeof c === "bigint") return [c, 1n];
  if (typeof c === "number") { if (!Number.isInteger(c)) throw new TypeError(`${c} is not a whole number; give fractions as Rationals or text such as "3/4"`); return [BigInt(c), 1n]; }
  if (typeof c === "string") { const [n, d] = c.replace("−", "-").split("/"); return normQ(BigInt(n.trim()), BigInt((d ?? "1").trim())); }
  if (c && typeof c.numerator === "bigint") return [c.numerator, c.denominator];
  throw new TypeError(`not a rational number: ${c}`);
}

function toIntegers(coeffs) {
  const pairs = coeffs.map(readQ);
  let L = 1n;
  for (const [, d] of pairs) L = lcmB(L, d);
  let ints = pairs.map(([n, d]) => n * (L / d));
  while (ints.length && ints[ints.length - 1] === 0n) ints.pop();
  if (ints.length < 2) throw new RangeError("the polynomial needs degree at least 1");
  let g = 0n;
  for (const c of ints) g = gcdB(g, c);
  ints = ints.map((c) => c / g);
  if (ints[ints.length - 1] < 0n) ints = ints.map((c) => -c);
  return ints;
}

const fmtPoly = (ints, v = "x") => {
  const parts = [];
  for (let i = ints.length - 1; i >= 0; i -= 1) {
    const c = ints[i];
    if (c === 0n) continue;
    const a = absB(c);
    const body = i === 0 ? `${a}` : `${a === 1n ? "" : a}${v}${i > 1 ? sup(i) : ""}`;
    parts.push((parts.length ? (c < 0n ? " − " : " + ") : c < 0n ? "−" : "") + body);
  }
  return parts.join("");
};

/**
 * Solves a polynomial over Q by radicals when that is possible.
 *
 * Returns { polynomial, degree, solvable, reason, factors } where each
 * factor is { polynomial, text, degree, multiplicity, group, order,
 * solvable, chain (the prime indices of the steps), roots }. For a solvable
 * factor `roots` lists { value (the exact root from @ratmath/algebraic),
 * expression (a RadicalExpression) } in the order real roots ascending, then
 * the others. Irreducible factors are handled up to degree 7.
 */
export function solveByRadicals(input, { expressions = true } = {}) {
  const f = toIntegers(readPolynomial(input));
  const { factors } = factorZ(f);
  const atoms = new Atoms();
  const out = [];
  for (const { factor: g, multiplicity } of factors) {
    const degree = g.length - 1;
    const text = fmtPoly(g);
    if (degree === 1) {
      const q = normQ(-g[0], g[1]);
      out.push({ polynomial: g, text, degree, multiplicity, group: "1", order: 1, solvable: true, chain: [], roots: [{ value: new Rational(q[0], q[1]), expression: new RadicalExpression(atoms, atoms.rat(q), null) }] });
      continue;
    }
    if (degree > 7) throw new RangeError(`${text} is irreducible of degree ${degree}; Galois groups are computed up to degree 7`);
    const action = galoisAction(g);
    const { ctx, rootBoxes } = polynomialContext(action);
    const entry = { polynomial: g, text, degree, multiplicity, group: action.name, order: action.order, solvable: Boolean(ctx), chain: ctx ? ctx.chain.map((s) => s.p) : null, roots: null };
    if (ctx && expressions) {
      entry.roots = withPrecision((P) => {
        ctx.setPrecision(P);
        const solver = new Solver(atoms);
        solver.P = P;
        const mark = atoms.mark();
        try {
          return action.roots.map((value, k) => {
            const expr = solver.express(ctx, ctx.gen(k), ctx.chain.length);
            const target = (F) => { const all = rootBoxes(F); return { box: all[k], others: all.filter((_, i) => i !== k) }; };
            return { value, expression: new RadicalExpression(atoms, expr, target) };
          });
        } catch (e) {
          if (e instanceof NeedPrecision) { atoms.rewind(mark); ctx.memo.clear(); }
          throw e;
        }
      });
    }
    out.push(entry);
  }
  const solvable = out.every((e) => e.solvable);
  const bad = out.find((e) => !e.solvable);
  const reason = solvable
    ? (out.length === 1 ? `the Galois group ${out[0].group} of order ${out[0].order} is solvable${out[0].chain.length ? `: it comes apart in steps of index ${out[0].chain.join(", ")}` : ""}` : "every irreducible factor has a solvable Galois group")
    : `the Galois group ${bad.group} of ${bad.text}, of order ${bad.order}, is not solvable: its chain of commutator subgroups stops before reaching the identity, so no tower of radicals reaches the roots`;
  return { polynomial: f, text: fmtPoly(f), degree: f.length - 1, solvable, reason, factors: out };
}

/** A primitive n-th root of unity e^(2πi/n) in radicals, for any n ≥ 1. */
export function rootOfUnityByRadicals(n) {
  if (!Number.isInteger(n) || n < 1) throw new RangeError("n is a positive whole number");
  const atoms = new Atoms();
  const expr = withPrecision((P) => {
    atoms.reset();
    const solver = new Solver(atoms);
    solver.P = P;
    // ζₙ = Π ζ_(p^e)^(cₚ) over the prime powers of n, with cₚ·(n/p^e) ≡ 1 (mod p^e)
    let out = atoms.rat(ONE);
    for (const p of new Set(primeFactorsSmall(n))) {
      let pe = 1;
      while (n % (pe * p) === 0) pe *= p;
      const m = n / pe;
      let c = 1;
      while ((c * m) % pe !== 1 % pe) c += 1;
      // ζ_(p^e): start from ζ_p and take p-th roots
      let z = solver.unity(p);
      for (let q = p * p; q <= pe; q *= p) z = solver.radical(p, z, [Math.cos(TWO_PI / q), Math.sin(TWO_PI / q)]);
      out = out.mul(z.pow(c));
    }
    return out;
  });
  const target = (F) => {
    const one = F.rational(ONE);
    const all = Array.from({ length: n }, (_, k) => F.rootNear(one, n, [Math.cos((TWO_PI * k) / n), Math.sin((TWO_PI * k) / n)]));
    return { box: all[1 % n], others: all.filter((_, k) => k !== 1 % n) };
  };
  return new RadicalExpression(atoms, expr, target);
}

/** cos(2π/n) in radicals, through the real subfield of the cyclotomic field; for a prime or 2·prime… any n whose odd part is prime or 1 when exact form is wanted, see regularPolygon. */
function cosByRadicals(p, atoms, solver) {
  // p an odd prime: (ζ + ζ⁻¹)/2 lies in the field fixed by ±1, the last step of the chain
  const ctx = cyclotomicContext(p);
  ctx.setPrecision(solver.P);
  const z = ctx.gen(0);
  const sum = ctx.add(z, ctx.act(p - 2, z)); // index p−2 is the element −1 ≡ p−1
  return solver.express(ctx, sum, ctx.chain.length).scale([1n, 2n]);
}

const totient = (n) => { let out = n; for (const p of new Set(primeFactorsSmall(n))) out = (out / p) * (p - 1); return out; };
const isPowerOfTwo = (n) => n >= 1 && (n & (n - 1)) === 0;

/**
 * The regular n-gon and ruler and compass: { n, constructible, totient,
 * reason, cos } where `cos` is cos(2π/n) as a RadicalExpression in real
 * square roots when the polygon is constructible (Gauss and Wantzel: exactly
 * when φ(n) is a power of 2), and in general radicals otherwise when n is
 * an odd prime.
 */
export function regularPolygon(n) {
  if (!Number.isInteger(n) || n < 3) throw new RangeError("a polygon has at least 3 sides");
  const phi = totient(n);
  const constructible = isPowerOfTwo(phi);
  const odd = [...new Set(primeFactorsSmall(n))].filter((p) => p !== 2);
  let twos = 0;
  for (let m = n; m % 2 === 0; m /= 2) twos += 1;
  const reason = constructible
    ? `φ(${n}) = ${phi} is a power of 2${!odd.length ? "" : !twos && odd.length === 1 ? `, ${n} being a Fermat prime` : `: ${n} = ${twos ? `2${twos > 1 ? sup(twos) : ""}·` : ""}${odd.join("·")} with ${odd.length > 1 ? "distinct Fermat primes" : "a Fermat prime"}`}, so cos(2π/${n}) is reached by square roots alone`
    : `φ(${n}) = ${phi} is not a power of 2, so cos(2π/${n}) has degree ${phi / 2} over Q and lies in no tower of square roots`;
  const oddPart = n >> twos;
  let cos = null;
  if (constructible || (odd.length === 1 && oddPart === odd[0] && twos <= 1)) {
    const atoms = new Atoms();
    const expr = withPrecision((P) => {
      atoms.reset();
      const solver = new Solver(atoms);
      solver.P = P;
      // cos and sin of 2π/m for the odd part m, by the addition formulas: 1/m ≡ Σ aₚ/p (mod 1)
      let c = atoms.rat(ONE);
      let s = new RPoly(atoms);
      for (const p of odd) {
        const rest = oddPart / p;
        let a = 1;
        while ((a * rest) % p !== 1 % p) a += 1;
        const cp = cosByRadicals(p, atoms, solver);
        const sp = solver.radical(2, atoms.rat(ONE).add(cp.mul(cp).neg()), [Math.sin(TWO_PI / p), 0]);
        for (let i = 0; i < a; i += 1) [c, s] = [c.mul(cp).add(s.mul(sp).neg()), s.mul(cp).add(c.mul(sp))];
      }
      // halve the angle once for each factor 2: cos(θ/2) = √((1 + cos θ)/2), with the sign of the true value
      let theta = TWO_PI / oddPart;
      for (let i = 0; i < twos; i += 1) {
        theta /= 2;
        c = solver.radical(2, c.add(atoms.rat(ONE)).scale([1n, 2n]), [Math.cos(theta), 0]);
      }
      return c;
    });
    const target = (F) => {
      // cos(2π/n) from the certified n-th root of unity
      const one = F.rational(ONE);
      const z = F.rootNear(one, n, [Math.cos(TWO_PI / n), Math.sin(TWO_PI / n)]);
      const zz = F.rootNear(one, n, [Math.cos((2 * TWO_PI) / n), Math.sin((2 * TWO_PI) / n)]);
      return { box: [z[0], z[1], 0n, 0n], others: n > 3 ? [[zz[0], zz[1], 0n, 0n]] : [] };
    };
    cos = new RadicalExpression(atoms, expr, target);
  }
  return { n, constructible, totient: phi, reason, cos };
}

/**
 * Can a root of the polynomial be constructed with ruler and compass from a
 * unit length? { constructible, reason, factors } with, per irreducible
 * factor, { text, degree, group, order, constructible, roots } and the roots
 * written in square roots when they can be. A number is constructible
 * exactly when the splitting field of its minimal polynomial has degree a
 * power of 2.
 */
export function constructible(input) {
  const f = toIntegers(readPolynomial(input));
  const { factors } = factorZ(f);
  const out = [];
  for (const { factor: g } of factors) {
    const degree = g.length - 1;
    const text = fmtPoly(g);
    if (!isPowerOfTwo(degree)) { out.push({ text, degree, group: null, order: null, constructible: false, reason: `its roots have degree ${degree} over Q, which is not a power of 2`, roots: null }); continue; }
    if (degree > 7) { out.push({ text, degree, group: null, order: null, constructible: null, reason: `degree ${degree} is a power of 2, but Galois groups are computed only up to degree 7, so this is left undecided`, roots: null }); continue; }
    const solved = solveByRadicals(g).factors[0];
    const ok = isPowerOfTwo(solved.order);
    out.push({
      text, degree, group: solved.group, order: solved.order, constructible: ok,
      reason: ok ? `its splitting field has degree ${solved.order}, a power of 2` : `its roots have degree ${degree}, a power of 2, but the splitting field has degree ${solved.order} (group ${solved.group}), which is not`,
      roots: ok ? solved.roots : null,
    });
  }
  const all = out.every((e) => e.constructible === true);
  const none = out.every((e) => e.constructible === false);
  return {
    polynomial: f,
    text: fmtPoly(f),
    constructible: all ? true : out.some((e) => e.constructible === null) ? null : out.some((e) => e.constructible) ? "some" : false,
    reason: out.length === 1 ? `${out[0].text}: ${out[0].reason}` : all ? "every root is constructible" : none ? "no root is constructible" : "some roots are constructible and some are not",
    factors: out,
  };
}

/** The three classical problems, each reduced to a polynomial and decided. */
export const classicalProblems = Object.freeze({
  /** Doubling the cube: constructing ∛2. */
  doublingTheCube() {
    const c = constructible("x^3 - 2");
    return { problem: "construct the side of a cube of twice the volume: a root of x³ − 2", possible: false, polynomial: "x³ − 2", reason: c.factors[0].reason };
  },
  /**
   * Trisecting an angle whose cosine is the rational number c (60° by
   * default): cos(θ/3) is a root of 4x³ − 3x − c.
   */
  trisecting(cosine = "1/2") {
    const q = readQ(cosine);
    if (absB(q[0]) > q[1]) throw new RangeError("a cosine lies between −1 and 1");
    // cos θ = 4c³ − 3c for c = cos(θ/3)
    const c = constructible([[-q[0], q[1]], -3n, 0n, 4n]);
    const cubic = c.factors.find((e) => e.degree === 3);
    return {
      problem: `trisect the angle whose cosine is ${q[0] < 0n ? "−" : ""}${fmtQ(q)}: construct a root of ${c.text}`,
      possible: !cubic,
      polynomial: c.text,
      reason: cubic ? `${cubic.text} is irreducible, so cos(θ/3) has degree 3 over Q, which is not a power of 2` : `${c.text} factors over Q into pieces of degree ${c.factors.map((e) => e.degree).join(" and ")}, so cos(θ/3) needs square roots at most`,
      roots: cubic ? null : c.factors.flatMap((e) => e.roots ?? []),
    };
  },
  /** Squaring the circle: constructing √π. Not a computation: π is transcendental (Lindemann, 1882). */
  squaringTheCircle() {
    return { problem: "construct a square with the area of a given circle: the length √π", possible: false, polynomial: null, reason: "every constructible number is algebraic, a root of a polynomial with rational coefficients, and π is not (Lindemann, 1882); this is a theorem quoted, not something computed here" };
  },
});
