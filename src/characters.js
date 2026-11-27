// Character tables of finite groups, exactly.
//
// The table is found by the Burnside–Dixon–Schneider method. The class sums
// of G span the centre of the group algebra, and each irreducible character χ
// gives a common eigenvector (|Cᵢ| χ(gᵢ) / χ(1))ᵢ of the matrices of class
// multiplication constants. Those eigenvectors are found over a prime field
// GF(p), with p ≡ 1 modulo the exponent e of G and p > 2√|G|, where every
// e-th root of unity exists and nothing has to be approximated. A value χ(g)
// is a sum of χ(1) roots of unity of order |g|; the number of times each one
// occurs is a small whole number recovered from the values mod p of χ on the
// powers of g (Dixon's formula). So the table comes out as exact elements of
// the cyclotomic field Q(ζ_e), and the orthogonality relations are then
// checked on the exact values.

import { cyclotomic } from "@ratmath/algebraic";
import { Group, factorInt, gcd } from "./group.js";

const CLASS_LIMIT = 400;
const SUB = "₀₁₂₃₄₅₆₇₈₉";
const SUP = "⁰¹²³⁴⁵⁶⁷⁸⁹";
const sub = (n) => String(n).replace(/\d/g, (d) => SUB[d]);
const sup = (n) => String(n).replace(/\d/g, (d) => SUP[d]);
const bgcd = (a, b) => { a = a < 0n ? -a : a; b = b < 0n ? -b : b; while (b) [a, b] = [b, a % b]; return a; };
const lcm = (a, b) => (a / gcd(a, b)) * b;
const mod = (a, n) => ((a % n) + n) % n;

// ---- cyclotomic numbers ------------------------------------------------------------------------

const FIELDS = new Map();
/** Q(ζₙ): the cyclotomic polynomial and every power of ζ reduced by it. */
function cycField(n) {
  let F = FIELDS.get(n);
  if (F) return F;
  const phi = cyclotomic(n).map(BigInt);
  const d = phi.length - 1;
  const pow = [Array.from({ length: d }, (_, i) => (i === 0 ? 1n : 0n))];
  for (let k = 1; k < n; k += 1) {
    const prev = pow[k - 1];
    const top = prev[d - 1];
    const next = Array(d).fill(0n);
    for (let i = 1; i < d; i += 1) next[i] = prev[i - 1];
    if (top !== 0n) for (let i = 0; i < d; i += 1) next[i] -= top * phi[i];
    pow.push(next);
  }
  const units = [];
  for (let k = 1; k <= n; k += 1) if (gcd(k, n) === 1) units.push(k % n);
  F = { n, degree: d, pow, units };
  FIELDS.set(n, F);
  return F;
}

/**
 * An element of a cyclotomic field Q(ζₙ), ζₙ = e^(2πi/n): the values of
 * characters. Held as whole-number coordinates in 1, ζ, …, ζ^(φ(n)−1) over a
 * common denominator, so equality is exact.
 */
export class Cyclotomic {
  constructor(n, c, den = 1n) {
    const F = cycField(n);
    let g = den;
    for (const v of c) if (v !== 0n) g = bgcd(g, v);
    if (g > 1n) { c = c.map((v) => v / g); den /= g; }
    this.n = n;
    this.c = c;
    this.den = den;
    this.F = F;
    Object.freeze(this);
  }
  static int(k, n = 1) { const F = cycField(n); const c = Array(F.degree).fill(0n); c[0] = BigInt(k); return new Cyclotomic(n, c); }
  /** ζₙᵏ. */
  static root(n, k = 1) { return new Cyclotomic(n, [...cycField(n).pow[mod(k, n)]]); }
  /** Σ mₖ ζₙᵏ from a list of multiplicities. */
  static fromPowers(n, m) {
    const F = cycField(n);
    const c = Array(F.degree).fill(0n);
    m.forEach((mk, k) => { if (mk) { const b = BigInt(mk); const p = F.pow[k % n]; for (let i = 0; i < F.degree; i += 1) c[i] += b * p[i]; } });
    return new Cyclotomic(n, c);
  }
  static of(x, n = 1) { return x instanceof Cyclotomic ? x : Cyclotomic.int(x, n); }

  /** The same number in the bigger field Q(ζₘ), n dividing m. */
  lift(m) {
    if (m === this.n) return this;
    if (m % this.n) throw new RangeError(`Q(ζ${sub(this.n)}) is not inside Q(ζ${sub(m)})`);
    const F = cycField(m);
    const step = m / this.n;
    const c = Array(F.degree).fill(0n);
    this.c.forEach((v, k) => { if (v !== 0n) { const p = F.pow[(k * step) % m]; for (let i = 0; i < F.degree; i += 1) c[i] += v * p[i]; } });
    return new Cyclotomic(m, c, this.den);
  }
  #with(o) {
    const b = Cyclotomic.of(o, this.n);
    if (b.n === this.n) return [this, b];
    const m = lcm(this.n, b.n);
    return [this.lift(m), b.lift(m)];
  }
  add(o) { const [a, b] = this.#with(o); return new Cyclotomic(a.n, a.c.map((v, i) => v * b.den + b.c[i] * a.den), a.den * b.den); }
  neg() { return new Cyclotomic(this.n, this.c.map((v) => -v), this.den); }
  sub(o) { return this.add(Cyclotomic.of(o, this.n).neg()); }
  mul(o) {
    const [a, b] = this.#with(o);
    const F = a.F;
    const c = Array(F.degree).fill(0n);
    for (let i = 0; i < F.degree; i += 1) {
      if (a.c[i] === 0n) continue;
      for (let j = 0; j < F.degree; j += 1) {
        if (b.c[j] === 0n) continue;
        const t = a.c[i] * b.c[j];
        const p = F.pow[(i + j) % F.n];
        for (let k = 0; k < F.degree; k += 1) c[k] += t * p[k];
      }
    }
    return new Cyclotomic(a.n, c, a.den * b.den);
  }
  /** Division by a whole number or fraction p/q given as (p, q). */
  over(p, q = 1) {
    let num = BigInt(q);
    let den = BigInt(p);
    if (den === 0n) throw new RangeError("division by zero");
    if (den < 0n) { den = -den; num = -num; }
    return new Cyclotomic(this.n, this.c.map((v) => v * num), this.den * den);
  }
  /** The Galois conjugate ζ ↦ ζᵏ (k coprime to n). */
  galois(k) {
    const F = this.F;
    const c = Array(F.degree).fill(0n);
    this.c.forEach((v, i) => { if (v !== 0n) { const p = F.pow[mod(i * k, F.n)]; for (let j = 0; j < F.degree; j += 1) c[j] += v * p[j]; } });
    return new Cyclotomic(this.n, c, this.den);
  }
  /** The complex conjugate. */
  conj() { return this.galois(-1); }
  equals(o) { const [a, b] = this.#with(o); return a.den === b.den && a.c.every((v, i) => v === b.c[i]); }
  get isZero() { return this.c.every((v) => v === 0n); }
  get isRational() { return this.c.every((v, i) => i === 0 || v === 0n); }
  get isInteger() { return this.isRational && this.den === 1n; }
  get isReal() { return this.equals(this.conj()); }
  /** A JavaScript number when the value is a whole number, else null. */
  toInteger() { return this.isInteger ? Number(this.c[0]) : null; }
  /** [numerator, denominator] as BigInts when rational, else null. */
  toFraction() { return this.isRational ? [this.c[0], this.den] : null; }
  /** { re, im } as floats. */
  approx() {
    let re = 0;
    let im = 0;
    const d = Number(this.den);
    this.c.forEach((v, k) => { if (v !== 0n) { const t = (2 * Math.PI * k) / this.n; re += Number(v) * Math.cos(t); im += Number(v) * Math.sin(t); } });
    return { re: re / d, im: im / d };
  }
  /** The distinct Galois conjugates of this number. */
  conjugates() {
    const out = [];
    for (const k of this.F.units) { const y = this.galois(k); if (!out.some((z) => z.equals(y))) out.push(y); }
    return out;
  }
  /** The smallest m with this number in Q(ζₘ). */
  get conductor() {
    const n = this.n;
    for (let m = 1; m <= n; m += 1) {
      if (n % m) continue;
      if (this.F.units.every((k) => k % m !== 1 % m || this.galois(k).equals(this))) return m;
    }
    return n;
  }
  /** Coordinates [numerators…], denominator in 1, ζₘ, …, for m the conductor. */
  #reduced() {
    const m = this.conductor;
    if (m === this.n) return { m, c: this.c, den: this.den };
    // solve Σ xⱼ ζₘʲ = this over Q, by elimination on whole numbers
    const big = this.F;
    const small = cycField(m);
    const step = this.n / m;
    const rows = Array.from({ length: big.degree }, (_, i) => [...Array.from({ length: small.degree }, (_, j) => big.pow[(j * step) % this.n][i]), this.c[i]]);
    const piv = [];
    let r = 0;
    for (let col = 0; col < small.degree; col += 1) {
      const p = rows.findIndex((row, i) => i >= r && row[col] !== 0n);
      if (p < 0) throw new Error("internal: conductor basis is dependent");
      [rows[r], rows[p]] = [rows[p], rows[r]];
      for (let i = 0; i < rows.length; i += 1) {
        if (i === r || rows[i][col] === 0n) continue;
        const a = rows[r][col];
        const b = rows[i][col];
        rows[i] = rows[i].map((v, k) => v * a - rows[r][k] * b);
        let g = 0n;
        for (const v of rows[i]) g = bgcd(g, v);
        if (g > 1n) rows[i] = rows[i].map((v) => v / g);
      }
      piv.push(r);
      r += 1;
    }
    // xⱼ = rhs/pivot·(1/den): put over one denominator
    let den = this.den;
    const fr = piv.map((ri, j) => [rows[ri][small.degree], rows[ri][j]]);
    let L = 1n;
    for (const [, q] of fr) { const a = q < 0n ? -q : q; L = (L / bgcd(L, a)) * a; }
    const c = fr.map(([p, q]) => (p * L) / q);
    den *= L;
    const out = new Cyclotomic(m, c, den);
    return { m, c: out.c, den: out.den };
  }
  #text(tex) {
    const frac = this.toFraction();
    const minus = tex ? "-" : "−";
    const showInt = (v) => (v < 0n ? `${minus}${-v}` : String(v));
    if (frac) return frac[1] === 1n ? showInt(frac[0]) : tex ? `${frac[0] < 0n ? "-" : ""}\\tfrac{${frac[0] < 0n ? -frac[0] : frac[0]}}{${frac[1]}}` : `${showInt(frac[0])}/${frac[1]}`;
    const conj = this.conjugates();
    if (conj.length === 2) {
      // a quadratic irrational: (t ± s√d)/q
      const t = this.add(conj[1]).toFraction();
      const nm = this.mul(conj[1]).toFraction();
      if (t && nm) {
        // 2x − t = ±√D with D = t² − 4·norm
        const q = t[1] * nm[1];
        const T = (t[0] * q) / t[1];
        const N = (nm[0] * q) / nm[1];
        let D = T * T - 4n * N * q; // (2q·x − T)² = D
        let s = 1n;
        const negative = D < 0n;
        if (negative) D = -D;
        for (let f = 2n; f * f <= D; f += 1n) while (D % (f * f) === 0n) { D /= f * f; s *= f; }
        // x = (T ± s√(±D)) / (2q)
        const a = this.approx();
        const signPlus = negative ? a.im > 0 : a.re * 2 * Number(q) - Number(T) > 0;
        let den = 2n * q;
        let g = bgcd(bgcd(T, s), den);
        const T2 = T / g;
        const s2 = s / g;
        den /= g;
        const rootText = negative ? (D === 1n ? "i" : tex ? `i\\sqrt{${D}}` : `i√${D}`) : tex ? `\\sqrt{${D}}` : `√${D}`;
        const sTerm = `${s2 === 1n ? "" : s2}${rootText}`;
        let body;
        if (T2 === 0n) body = `${signPlus ? "" : minus}${sTerm}`;
        else body = `${showInt(T2)} ${signPlus ? "+" : minus} ${sTerm}`;
        if (den === 1n) return body;
        if (tex) return T2 === 0n ? `${signPlus ? "" : "-"}\\tfrac{${sTerm}}{${den}}` : `\\tfrac{${body}}{${den}}`;
        return T2 === 0n ? `${body}/${den}` : `(${body})/${den}`;
      }
    }
    const red = this.#reduced();
    const { m, den } = red;
    let c = red.c;
    if (m > 2 && factorInt(m).size === 1 && factorInt(m).get(m) === 1) {
      // prime conductor: 1 = −(ζ + … + ζ^(m−1)), so the constant can be traded for the last power
      const alt = [0n, ...c.slice(1).map((v) => v - c[0]), -c[0]];
      if (alt.filter((v) => v !== 0n).length < c.filter((v) => v !== 0n).length) c = alt;
    } else if (den === 1n) {
      // a single root of unity, up to sign
      for (let k = 1; k < m; k += 1) {
        for (const sgn of [1n, -1n]) {
          if (cycField(m).pow[k].every((v, i) => v * sgn === c[i])) { c = Array(k + 1).fill(0n); c[k] = sgn; }
        }
        if (c.length === k + 1 && c[k] !== 0n && c.slice(0, k).every((v) => v === 0n)) break;
      }
    }
    const z = tex ? `\\zeta_{${m}}` : `ζ${sub(m)}`;
    let s = "";
    c.forEach((v, k) => {
      if (v === 0n) return;
      const a = v < 0n ? -v : v;
      const mon = k === 0 ? String(a) : `${a === 1n ? "" : a}${z}${k === 1 ? "" : tex ? `^{${k}}` : sup(k)}`;
      s += s ? ` ${v < 0n ? minus : "+"} ${mon}` : `${v < 0n ? minus : ""}${mon}`;
    });
    if (den === 1n) return s;
    return tex ? `\\tfrac{${s}}{${den}}` : `(${s})/${den}`;
  }
  toString() { return this.#text(false); }
  toTeX() { return this.#text(true); }
}

// ---- arithmetic mod p ----------------------------------------------------------------------------

function powMod(a, e, p) {
  let r = 1;
  a %= p;
  while (e > 0) {
    if (e & 1) r = (r * a) % p;
    a = (a * a) % p;
    e = Math.floor(e / 2);
  }
  return r;
}
const invMod = (a, p) => powMod(mod(a, p), p - 2, p);
const isPrimeN = (n) => { if (n < 2) return false; for (let d = 2; d * d <= n; d += 1) if (n % d === 0) return false; return true; };
function primitiveRootMod(p) {
  const qs = [...factorInt(p - 1).keys()];
  for (let g = 2; g < p; g += 1) if (qs.every((q) => powMod(g, (p - 1) / q, p) !== 1)) return g;
  return 1;
}

/** Row-reduce vectors mod p: { basis (reduced echelon), pivots }. */
function rref(rows, p) {
  const basis = [];
  const pivots = [];
  const width = rows.length ? rows[0].length : 0;
  const work = rows.map((r) => Float64Array.from(r));
  let r = 0;
  for (let col = 0; col < width && r < work.length; col += 1) {
    let piv = -1;
    for (let i = r; i < work.length; i += 1) if (work[i][col] !== 0) { piv = i; break; }
    if (piv < 0) continue;
    [work[r], work[piv]] = [work[piv], work[r]];
    const inv = invMod(work[r][col], p);
    const row = work[r];
    for (let k = col; k < width; k += 1) row[k] = (row[k] * inv) % p;
    for (let i = 0; i < work.length; i += 1) {
      if (i === r) continue;
      const f = work[i][col];
      if (f === 0) continue;
      const wi = work[i];
      for (let k = col; k < width; k += 1) if (row[k] !== 0) wi[k] = mod(wi[k] - f * row[k], p);
    }
    pivots.push(col);
    r += 1;
  }
  for (let i = 0; i < r; i += 1) basis.push(work[i]);
  return { basis, pivots };
}

/** The null space of a square matrix (rows) mod p, as column coordinate vectors. */
function nullSpaceMod(A, p) {
  const m = A.length;
  const { basis, pivots } = rref(A, p);
  const free = [];
  for (let c = 0; c < m; c += 1) if (!pivots.includes(c)) free.push(c);
  return free.map((f) => {
    const v = new Float64Array(m);
    v[f] = 1;
    pivots.forEach((pc, i) => { v[pc] = mod(-basis[i][f], p); });
    return v;
  });
}

/** The characteristic polynomial of a matrix mod p (constant first), by Hessenberg form. */
function charPolyMod(Ain, p) {
  const m = Ain.length;
  const H = Ain.map((r) => Float64Array.from(r));
  for (let k = 1; k < m - 1; k += 1) {
    let piv = -1;
    for (let i = k; i < m; i += 1) if (H[i][k - 1] !== 0) { piv = i; break; }
    if (piv < 0) continue;
    if (piv !== k) {
      [H[k], H[piv]] = [H[piv], H[k]];
      for (let i = 0; i < m; i += 1) [H[i][k], H[i][piv]] = [H[i][piv], H[i][k]];
    }
    const inv = invMod(H[k][k - 1], p);
    for (let i = k + 1; i < m; i += 1) {
      const f = (H[i][k - 1] * inv) % p;
      if (f === 0) continue;
      for (let j = 0; j < m; j += 1) H[i][j] = mod(H[i][j] - f * H[k][j], p);
      for (let j = 0; j < m; j += 1) H[j][k] = (H[j][k] + f * H[j][i]) % p;
    }
  }
  // pₖ(x) = (x − hₖₖ) pₖ₋₁ − Σᵢ hₖ₋ᵢ,ₖ (Π subdiagonals) pₖ₋ᵢ₋₁
  const polys = [[1]];
  for (let k = 1; k <= m; k += 1) {
    const prev = polys[k - 1];
    const cur = new Array(k + 1).fill(0);
    for (let i = 0; i < prev.length; i += 1) {
      cur[i + 1] = (cur[i + 1] + prev[i]) % p;
      cur[i] = mod(cur[i] - H[k - 1][k - 1] * prev[i], p);
    }
    let prod = 1;
    for (let i = 1; i < k; i += 1) {
      prod = (prod * H[k - i][k - i - 1]) % p;
      if (prod === 0) break;
      const f = (H[k - i - 1][k - 1] * prod) % p;
      if (f === 0) continue;
      const q = polys[k - i - 1];
      for (let j = 0; j < q.length; j += 1) cur[j] = mod(cur[j] - f * q[j], p);
    }
    polys.push(cur);
  }
  return polys[m];
}

function rootsMod(f, p) {
  const out = [];
  const d = f.length - 1;
  for (let x = 0; x < p && out.length < d; x += 1) {
    let v = 0;
    for (let i = d; i >= 0; i -= 1) v = (v * x + f[i]) % p;
    if (v === 0) out.push(x);
  }
  return out;
}

// ---- class functions -----------------------------------------------------------------------------

/**
 * A class function on a group: one cyclotomic value per conjugacy class, in
 * the order of the table's classes. Characters are the class functions that
 * are whole-number combinations of the irreducible ones with no negative
 * multiplicity.
 */
export class ClassFunction {
  constructor(table, values, name = null) {
    if (values.length !== table.classes.length) throw new RangeError(`expected ${table.classes.length} values, one per class`);
    this.table = table;
    this.values = Object.freeze(values.map((v) => Cyclotomic.of(v)));
    this.name = name;
    Object.freeze(this);
  }
  /** The value at the identity: a number when whole, else the cyclotomic value. */
  get degree() { const v = this.values[0]; return v.toInteger() ?? v; }
  value(x) { return this.values[this.table.classOf(x)]; }
  #same(o) { if (!(o instanceof ClassFunction) || o.table !== this.table) throw new TypeError("expected a class function on the same table"); return o; }
  add(o) { this.#same(o); return new ClassFunction(this.table, this.values.map((v, i) => v.add(o.values[i]))); }
  sub(o) { this.#same(o); return new ClassFunction(this.table, this.values.map((v, i) => v.sub(o.values[i]))); }
  /** The product: the character of the tensor product. */
  mul(o) { this.#same(o); return new ClassFunction(this.table, this.values.map((v, i) => v.mul(o.values[i]))); }
  tensor(o) { return this.mul(o); }
  scale(k) { return new ClassFunction(this.table, this.values.map((v) => v.mul(k))); }
  pow(k) { let r = this.table.trivial; for (let i = 0; i < k; i += 1) r = r.mul(this); return r; }
  /** The complex conjugate: the character of the dual representation. */
  conj() { return new ClassFunction(this.table, this.values.map((v) => v.conj())); }
  dual() { return this.conj(); }
  equals(o) { this.#same(o); return this.values.every((v, i) => v.equals(o.values[i])); }
  /** ⟨χ, ψ⟩ = (1/|G|) Σ χ(g) ψ(g)*: a number when whole, else a Cyclotomic. */
  inner(o) {
    this.#same(o);
    let s = Cyclotomic.int(0);
    this.table.classes.forEach((c, i) => { s = s.add(this.values[i].mul(o.values[i].conj()).mul(c.size)); });
    const v = s.over(this.table.group.order);
    return v.toInteger() ?? v;
  }
  get norm() { return this.inner(this); }
  /** The multiplicity of each irreducible character, in table order. */
  multiplicities() { return this.table.irreducibles.map((chi) => this.inner(chi)); }
  get isCharacter() { return this.multiplicities().every((m) => Number.isInteger(m) && m >= 0) && !this.values.every((v) => v.isZero); }
  get isVirtualCharacter() { return this.multiplicities().every((m) => Number.isInteger(m)); }
  get isIrreducible() { return this.norm === 1 && this.table.irreducibles.some((chi) => chi.equals(this)); }
  get isLinear() { return this.degree === 1 && this.isIrreducible; }
  /** [{ character, multiplicity }] for the irreducible constituents. */
  decompose() {
    const m = this.multiplicities();
    return this.table.irreducibles.map((character, i) => ({ character, multiplicity: m[i] })).filter((t) => t.multiplicity !== 0);
  }
  /** "χ₁ + 2χ₃": the decomposition as text. */
  decomposition() {
    const parts = this.decompose().map(({ character, multiplicity }) => `${multiplicity === 1 ? "" : multiplicity === -1 ? "−" : String(multiplicity).replace("-", "−")}${character.name}`);
    return parts.length ? parts.join(" + ").replace(/\+ −/g, "− ") : "0";
  }
  /** g ↦ χ(gᵏ): the k-th Adams operation. */
  adams(k) { const pm = this.table.powerMap(k); return new ClassFunction(this.table, pm.map((j) => this.values[j])); }
  /** The character of the symmetric square: (χ(g)² + χ(g²))/2. */
  symmetricSquare() { const a = this.adams(2); return new ClassFunction(this.table, this.values.map((v, i) => v.mul(v).add(a.values[i]).over(2))); }
  /** The character of the alternating square: (χ(g)² − χ(g²))/2. */
  alternatingSquare() { const a = this.adams(2); return new ClassFunction(this.table, this.values.map((v, i) => v.mul(v).sub(a.values[i]).over(2))); }
  /** The Frobenius–Schur indicator (1/|G|) Σ χ(g²): 1, 0 or −1 for an irreducible character. */
  indicator() { return this.adams(2).inner(this.table.trivial); }
  /** "real", "complex" or "quaternionic" for an irreducible character, by its indicator. */
  get type() { const v = this.indicator(); return v === 1 ? "real" : v === 0 ? "complex" : v === -1 ? "quaternionic" : null; }
  /** Indices of the classes where χ(g) = χ(1). */
  kernelClasses() { return this.values.map((v, i) => (v.equals(this.values[0]) ? i : -1)).filter((i) => i >= 0); }
  /** The kernel {g : χ(g) = χ(1)}, a normal subgroup of the group. */
  kernel() { return this.table.subgroupOfClasses(this.kernelClasses()); }
  get isFaithful() { return this.kernelClasses().length === 1; }
  /** Classes where |χ(g)| = χ(1): the centre Z(χ), the elements acting as scalars. */
  centreClasses() { const d2 = this.values[0].mul(this.values[0]); return this.values.map((v, i) => (v.mul(v.conj()).equals(d2) ? i : -1)).filter((i) => i >= 0); }
  /** The restriction to a subgroup H, as a class function on H's table. */
  restrict(H) {
    const T = H instanceof CharacterTable ? H : characterTable(H);
    return new ClassFunction(T, T.classes.map((c) => this.value(c.representative)));
  }
  /** The character induced from this one (on a subgroup) up to the group of `table`. */
  induce(table) {
    const T = table instanceof CharacterTable ? table : characterTable(table);
    const H = this.table;
    const G = T.group;
    const sums = T.classes.map(() => Cyclotomic.int(0));
    H.classes.forEach((c, i) => {
      const k = T.classOf(c.representative);
      sums[k] = sums[k].add(this.values[i].over(H.group.order / c.size));
    });
    return new ClassFunction(T, sums.map((s, k) => s.mul(G.order / T.classes[k].size)));
  }
  toString() { return `${this.name ? `${this.name} = ` : ""}[${this.values.join(", ")}]`; }
}

// ---- the table -----------------------------------------------------------------------------------

const TABLES = new WeakMap();

/** The character table of a finite group (cached on the group). */
export function characterTable(G, options = {}) {
  if (!(G instanceof Group)) throw new TypeError("characterTable needs a finite group");
  let T = TABLES.get(G);
  if (!T) { T = new CharacterTable(G, options); TABLES.set(G, T); }
  return T;
}

export class CharacterTable {
  #classOfIndex;
  #powerMaps = new Map();
  #irreducibles;

  constructor(G, { classLimit = CLASS_LIMIT } = {}) {
    this.group = G;
    const cl = G.conjugacyClasses();
    if (cl.length > classLimit) throw new RangeError(`${cl.length} conjugacy classes: character tables are limited to ${classLimit}`);
    const letters = new Map();
    this.classes = Object.freeze(cl.map((c) => {
      const k = letters.get(c.order) ?? 0;
      letters.set(c.order, k + 1);
      return Object.freeze({ name: `${c.order}${k < 26 ? String.fromCharCode(97 + k) : `_${k + 1}`}`, representative: c.representative, index: c.indices[0], size: c.size, order: c.order, centralizerOrder: G.order / c.size });
    }));
    this.#classOfIndex = new Int32Array(G.order);
    cl.forEach((c, i) => { for (const x of c.indices) this.#classOfIndex[x] = i; });
    this.exponent = cl.reduce((e, c) => lcm(e, c.order), 1);
    this.inverseClass = Object.freeze(this.classes.map((c) => this.#classOfIndex[G.invIndex(c.index)]));
    this.#compute();
    Object.freeze(this);
  }

  /** The index of the class of an element. */
  classOf(x) {
    const i = this.group.indexOf(x);
    if (i < 0) throw new RangeError("not an element of the group");
    return this.#classOfIndex[i];
  }
  /** For each class, the class of the k-th powers of its elements. */
  powerMap(k) {
    const e = this.exponent;
    const kk = mod(k, e);
    let pm = this.#powerMaps.get(kk);
    if (!pm) {
      pm = Object.freeze(this.classes.map((c) => this.#classOfIndex[this.group.powIndex(c.index, mod(k, c.order))]));
      this.#powerMaps.set(kk, pm);
    }
    return pm;
  }

  #compute() {
    const G = this.group;
    const r = this.classes.length;
    const e = this.exponent;
    const cls = this.#classOfIndex;
    const sizes = this.classes.map((c) => c.size);
    let p = e + 1;
    while (p * p <= 4 * G.order || !isPrimeN(p)) p += e;
    this.prime = p;

    // class matrix Mⱼ[i][k] = #{y ∈ Cⱼ : gₖ y⁻¹ ∈ Cᵢ}, so that Cᵢ·Cⱼ = Σₖ Mⱼ[i][k] Cₖ
    const members = this.classes.map(() => []);
    for (let x = 0; x < G.order; x += 1) members[cls[x]].push(x);
    const classMatrix = (j) => {
      const M = Array.from({ length: r }, () => new Float64Array(r));
      for (let k = 0; k < r; k += 1) {
        const gk = this.classes[k].index;
        for (const y of members[j]) M[cls[G.mulIndex(gk, G.invIndex(y))]][k] += 1;
      }
      for (const row of M) for (let k = 0; k < r; k += 1) row[k] %= p;
      return M;
    };

    // split the whole space into common eigenspaces
    let spaces = [rref(Array.from({ length: r }, (_, i) => { const v = new Float64Array(r); v[i] = 1; return v; }), p)];
    const order = [...Array(r).keys()].slice(1).sort((a, b) => sizes[a] - sizes[b]);
    this.matricesUsed = 0;
    for (const j of order) {
      if (spaces.every((S) => S.basis.length === 1)) break;
      const M = classMatrix(j);
      this.matricesUsed += 1;
      const next = [];
      for (const S of spaces) {
        const m = S.basis.length;
        if (m === 1) { next.push(S); continue; }
        // images of the basis vectors, and the matrix of Mⱼ in that basis (read off at the pivots)
        const A = Array.from({ length: m }, () => new Float64Array(m));
        for (let c = 0; c < m; c += 1) {
          const b = S.basis[c];
          for (let i = 0; i < m; i += 1) {
            const row = M[S.pivots[i]];
            let s = 0;
            for (let k = 0; k < r; k += 1) if (b[k] !== 0 && row[k] !== 0) s = (s + row[k] * b[k]) % p;
            A[i][c] = s;
          }
        }
        const roots = rootsMod(charPolyMod(A, p), p);
        if (roots.length === 1) { next.push(S); continue; }
        let total = 0;
        for (const lam of roots) {
          const B = A.map((row, i) => { const t = Float64Array.from(row); t[i] = mod(t[i] - lam, p); return t; });
          const ns = nullSpaceMod(B, p);
          total += ns.length;
          const vecs = ns.map((co) => {
            const v = new Float64Array(r);
            for (let c = 0; c < m; c += 1) if (co[c] !== 0) { const b = S.basis[c]; for (let k = 0; k < r; k += 1) if (b[k] !== 0) v[k] = (v[k] + co[c] * b[k]) % p; }
            return v;
          });
          next.push(rref(vecs, p));
        }
        if (total !== m) throw new Error("internal: a class matrix did not split into eigenspaces");
      }
      spaces = next;
    }
    if (spaces.length !== r || spaces.some((S) => S.basis.length !== 1)) throw new Error("internal: the class matrices did not separate the characters");

    // each line gives a character mod p: ωᵢ = |Cᵢ| χ(gᵢ)/χ(1), normalised by ω₁ = 1
    const w = primitiveRootMod(p);
    const omega = powMod(w, (p - 1) / e, p);
    const maxDeg = Math.floor(Math.sqrt(G.order));
    const chars = spaces.map((S) => {
      const v = S.basis[0];
      const s0 = invMod(v[0], p);
      const om = Array.from(v, (x) => (x * s0) % p);
      let s = 0;
      for (let i = 0; i < r; i += 1) s = (s + ((om[i] * om[this.inverseClass[i]]) % p) * invMod(sizes[i], p)) % p;
      // s = |G|/d²
      const d2 = (mod(G.order, p) * invMod(s, p)) % p;
      let d = 0;
      for (let t = 1; t <= maxDeg; t += 1) if (G.order % t === 0 && (t * t) % p === d2) { d = t; break; }
      if (!d) throw new Error("internal: no degree found for a character");
      const modp = om.map((x, i) => (((x * d) % p) * invMod(sizes[i], p)) % p);
      // Dixon: the multiplicity of εᵏ among the eigenvalues of g is (1/n) Σₗ χ(gˡ) ε^(−kl)
      const values = this.classes.map((c, i) => {
        const n = c.order;
        if (n === 1) return Cyclotomic.int(d, e);
        const eps = powMod(omega, e / n, p);
        const epsInv = invMod(eps, p);
        const along = Array.from({ length: n }, (_, l) => modp[cls[G.powIndex(c.index, l)]]);
        const ninv = invMod(n, p);
        const mult = new Array(e).fill(0);
        let count = 0;
        for (let k = 0; k < n; k += 1) {
          const z = powMod(epsInv, k, p);
          let acc = 0;
          let zl = 1;
          for (let l = 0; l < n; l += 1) { acc = (acc + along[l] * zl) % p; zl = (zl * z) % p; }
          const mk = (acc * ninv) % p;
          if (mk > d) throw new Error("internal: a multiplicity came out too large");
          mult[k * (e / n)] = mk;
          count += mk;
        }
        if (count !== d) throw new Error("internal: multiplicities do not add up to the degree");
        return Cyclotomic.fromPowers(e, mult);
      });
      return { d, values };
    });

    // order: by degree, the trivial character first, then real values before complex, larger before smaller
    const keyOf = (ch) => ch.values.map((v) => { const a = v.approx(); return [a.re, a.im]; });
    chars.forEach((ch) => { ch.key = keyOf(ch); ch.trivial = ch.values.every((v) => v.toInteger() === 1); });
    chars.sort((a, b) => {
      if (a.d !== b.d) return a.d - b.d;
      if (a.trivial !== b.trivial) return a.trivial ? -1 : 1;
      for (let i = 0; i < r; i += 1) {
        if (Math.abs(a.key[i][0] - b.key[i][0]) > 1e-9) return b.key[i][0] - a.key[i][0];
        if (Math.abs(a.key[i][1] - b.key[i][1]) > 1e-9) return b.key[i][1] - a.key[i][1];
      }
      return 0;
    });
    this.#irreducibles = Object.freeze(chars.map((ch, i) => new ClassFunction(this, ch.values, `χ${sub(i + 1)}`)));
    const check = this.verify();
    if (!check.ok) throw new Error("internal: the character table failed its orthogonality check");
  }

  /** The irreducible characters, by degree, the trivial one first. */
  get irreducibles() { return this.#irreducibles; }
  get characters() { return this.#irreducibles; }
  get degrees() { return this.#irreducibles.map((c) => c.degree); }
  get trivial() { return this.#irreducibles[0]; }
  get linearCharacters() { return this.#irreducibles.filter((c) => c.degree === 1); }
  /** The regular character: |G| at the identity and 0 elsewhere. */
  get regular() { return new ClassFunction(this, this.classes.map((c, i) => (i === 0 ? this.group.order : 0))); }
  /** A class function from a list of values (numbers or Cyclotomics) or a function of a representative. */
  classFunction(values) {
    if (typeof values === "function") return new ClassFunction(this, this.classes.map((c) => values(c.representative, c)));
    return new ClassFunction(this, values);
  }
  /** The character of a permutation action (the number of fixed points) or of the cosets of a subgroup. */
  permutationCharacter(X) {
    if (X instanceof Group) return characterTable(X).trivial.induce(this);
    if (X && typeof X.fixedPoints === "function") return new ClassFunction(this, this.classes.map((c) => { const f = X.fixedPoints(c.representative); return typeof f === "number" ? f : f.length; }));
    throw new TypeError("expected an action of the group or a subgroup");
  }
  /** Decompose a class function: the same as cf.decompose(). */
  decompose(cf) { return cf.decompose(); }

  /**
   * Check the table on its exact values: row orthogonality ⟨χᵢ, χⱼ⟩ = δᵢⱼ,
   * column orthogonality Σ χ(g)χ(h)* = |C(g)| or 0, and Σ χ(1)² = |G|.
   */
  verify() {
    const irr = this.#irreducibles;
    const r = irr.length;
    let rows = true;
    for (let i = 0; i < r && rows; i += 1) for (let j = i; j < r; j += 1) if (irr[i].inner(irr[j]) !== (i === j ? 1 : 0)) { rows = false; break; }
    let columns = true;
    for (let a = 0; a < r && columns; a += 1) {
      for (let b = a; b < r; b += 1) {
        let s = Cyclotomic.int(0);
        for (const chi of irr) s = s.add(chi.values[a].mul(chi.values[b].conj()));
        if (s.toInteger() !== (a === b ? this.classes[a].centralizerOrder : 0)) { columns = false; break; }
      }
    }
    const degrees = irr.reduce((s, c) => s + c.degree * c.degree, 0) === this.group.order && irr.every((c) => this.group.order % c.degree === 0);
    return { rows, columns, degrees, ok: rows && columns && degrees };
  }

  /** The subgroup made of some whole classes (given by index). */
  subgroupOfClasses(list) {
    const G = this.group;
    const set = new Set(list);
    const indices = [];
    for (let x = 0; x < G.order; x += 1) if (set.has(this.#classOfIndex[x])) indices.push(x);
    return G.subgroupFromIndices(indices);
  }
  /**
   * The normal subgroups read off the table: every one is an intersection of
   * kernels of irreducible characters. [{ order, classes, kernelOf }], by order.
   */
  normalSubgroups() {
    const r = this.classes.length;
    const kernels = this.#irreducibles.map((c) => c.kernelClasses());
    const seen = new Map();
    const all = [...Array(r).keys()];
    const add = (set, from) => {
      const key = set.join(",");
      if (seen.has(key)) return false;
      seen.set(key, { classes: set, order: set.reduce((s, i) => s + this.classes[i].size, 0), kernelOf: from });
      return true;
    };
    add(all, [0]);
    let frontier = [all];
    while (frontier.length) {
      const next = [];
      for (const S of frontier) {
        kernels.forEach((K, i) => {
          const T = S.filter((c) => K.includes(c));
          if (T.length < S.length && add(T, [...new Set([...seen.get(S.join(",")).kernelOf, i])])) next.push(T);
        });
      }
      frontier = next;
    }
    return [...seen.values()].sort((a, b) => a.order - b.order).map((n) => Object.freeze({ order: n.order, classes: n.classes.map((i) => this.classes[i].name), kernelOf: n.kernelOf.map((i) => this.#irreducibles[i].name) }));
  }
  get isSimple() { return this.group.order > 1 && this.normalSubgroups().length === 2; }
  get isAbelian() { return this.degrees.every((d) => d === 1); }
  /** |G/G′|, the number of linear characters. */
  get abelianizationOrder() { return this.linearCharacters.length; }
  /** |G′|: the commutator subgroup is the intersection of the kernels of the linear characters. */
  get derivedSubgroupOrder() { return this.group.order / this.abelianizationOrder; }
  /** |Z(G)|: the number of classes of size 1. */
  get centreOrder() { return this.classes.filter((c) => c.size === 1).length; }
  /** Classes equal to their inverse class; their number is the number of real-valued irreducible characters. */
  get realClasses() { return this.classes.filter((c, i) => this.inverseClass[i] === i).map((c) => c.name); }

  /** Rows of text: the header and one row per character. */
  rows({ tex = false } = {}) {
    const head = ["", ...this.classes.map((c) => c.name)];
    const sizes = ["size", ...this.classes.map((c) => String(c.size))];
    const body = this.#irreducibles.map((chi) => [tex ? `\\chi_{${chi.name.slice(1).replace(/[₀-₉]/g, (d) => SUB.indexOf(d))}}` : chi.name, ...chi.values.map((v) => (tex ? v.toTeX() : String(v)))]);
    return [head, sizes, ...body];
  }
  toString() {
    const rows = this.rows();
    const width = rows[0].map((_, j) => Math.max(...rows.map((row) => [...row[j]].length)));
    return rows.map((row) => row.map((cell, j) => (j === 0 ? cell.padEnd(width[0]) : " ".repeat(width[j] - [...cell].length) + cell)).join("  ")).join("\n");
  }
  toTeX() {
    const rows = this.rows({ tex: true });
    return `\\begin{array}{c|${"r".repeat(this.classes.length)}}\n${rows.map((row, i) => row.join(" & ") + (i === 1 ? " \\\\ \\hline" : " \\\\")).join("\n")}\n\\end{array}`;
  }
  describe() {
    const G = this.group;
    return `${G.name ?? "group"} of order ${G.order}: ${this.classes.length} classes, degrees ${this.degrees.join(", ")}; values in Q(ζ${sub(this.exponent)}), found mod ${this.prime}`;
  }
}
