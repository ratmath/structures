// Towers of field extensions of Q.
//
// extension({ a: "a^2 - 2", b: "b^3 - a" }) is Q(a, b) with a² = 2 and
// b³ = a: each polynomial may use the generators before it and must be
// irreducible over the field they generate, which is checked. Elements are
// written in the tower's own basis (the products aⁱbʲ with i < 2, j < 3),
// and the degree is the product of the steps: the tower law.
//
// Underneath, the whole tower is one number field Q(θ) from
// @ratmath/algebraic (a primitive element found at each step by Trager's
// method), which is where products, inverses and factoring happen. The
// change of basis between the powers of θ and the tower basis is an exact
// matrix, so nothing is lost either way.

import { Rational } from "@ratmath/core";
import { NumberField, factorOver, extend, rootsIn, roots as algebraicRoots } from "@ratmath/algebraic";
import { ring } from "./mpoly.js";
import { fractions } from "./ratfun.js";
import { Group } from "./group.js";

// ---- rational numbers as BigInt pairs, and the little linear algebra needed

const gcdB = (a, b) => { a = a < 0n ? -a : a; b = b < 0n ? -b : b; while (b) [a, b] = [b, a % b]; return a; };
const norm = (n, d) => { if (d < 0n) { n = -n; d = -d; } const g = gcdB(n, d) || 1n; return [n / g, d / g]; };
const Z = [0n, 1n];
const ONE = [1n, 1n];
const addQ = (a, b) => (a[0] === 0n ? b : b[0] === 0n ? a : norm(a[0] * b[1] + b[0] * a[1], a[1] * b[1]));
const subQ = (a, b) => addQ(a, [-b[0], b[1]]);
const mulQ = (a, b) => (a[0] === 0n || b[0] === 0n ? Z : norm(a[0] * b[0], a[1] * b[1]));
const divQ = (a, b) => norm(a[0] * b[1], a[1] * b[0]);
const toQ = (x) => {
  if (Array.isArray(x)) return norm(BigInt(x[0]), BigInt(x[1]));
  if (typeof x === "bigint") return [x, 1n];
  if (typeof x === "number") { if (!Number.isInteger(x)) throw new TypeError(`${x} is not a whole number; write fractions as text such as "3/4"`); return [BigInt(x), 1n]; }
  if (x instanceof Rational) return [x.numerator, x.denominator];
  if (typeof x === "string" && /^\s*-?\d+\s*(\/\s*\d+\s*)?$/.test(x)) { const [n, d] = x.split("/"); return norm(BigInt(n.trim()), BigInt((d ?? "1").trim())); }
  if (x && typeof x.numerator === "bigint") return [x.numerator, x.denominator];
  throw new TypeError(`not a rational number: ${x}`);
};
const pub = ([n, d]) => new Rational(n, d);

/** Solve A·x = b over Q (A as a list of columns); null when there is no solution, the first solution otherwise. */
function solveColumns(cols, b) {
  const rows = b.length;
  const n = cols.length;
  const M = Array.from({ length: rows }, (_, i) => [...cols.map((c) => c[i]), b[i]]);
  const pivots = [];
  let r = 0;
  for (let c = 0; c < n && r < rows; c += 1) {
    let p = r;
    while (p < rows && M[p][c][0] === 0n) p += 1;
    if (p === rows) continue;
    [M[r], M[p]] = [M[p], M[r]];
    const inv = divQ(ONE, M[r][c]);
    M[r] = M[r].map((x) => mulQ(x, inv));
    for (let i = 0; i < rows; i += 1) {
      if (i === r || M[i][c][0] === 0n) continue;
      const f = M[i][c];
      M[i] = M[i].map((x, j) => subQ(x, mulQ(f, M[r][j])));
    }
    pivots.push(c);
    r += 1;
  }
  for (let i = r; i < rows; i += 1) if (M[i][n][0] !== 0n) return null;
  const x = Array(n).fill(Z);
  pivots.forEach((c, i) => { x[c] = M[i][n]; });
  return x;
}

/** A basis of the null space of the matrix with the given rows. */
function nullSpace(rowsIn, n) {
  const M = rowsIn.map((r) => [...r]);
  const pivots = [];
  let r = 0;
  for (let c = 0; c < n && r < M.length; c += 1) {
    let p = r;
    while (p < M.length && M[p][c][0] === 0n) p += 1;
    if (p === M.length) continue;
    [M[r], M[p]] = [M[p], M[r]];
    const inv = divQ(ONE, M[r][c]);
    M[r] = M[r].map((x) => mulQ(x, inv));
    for (let i = 0; i < M.length; i += 1) {
      if (i === r || M[i][c][0] === 0n) continue;
      const f = M[i][c];
      M[i] = M[i].map((x, j) => subQ(x, mulQ(f, M[r][j])));
    }
    pivots.push(c);
    r += 1;
  }
  const free = [];
  for (let c = 0; c < n; c += 1) if (!pivots.includes(c)) free.push(c);
  return free.map((f) => {
    const v = Array(n).fill(Z);
    v[f] = ONE;
    pivots.forEach((c, i) => { v[c] = [-M[i][f][0], M[i][f][1]]; });
    return v;
  });
}

const SUP = "⁰¹²³⁴⁵⁶⁷⁸⁹";
const sup = (k) => String(k).split("").map((c) => SUP[Number(c)]).join("");
const fmtQ = ([n, d]) => (d === 1n ? `${n < 0n ? -n : n}` : `${n < 0n ? -n : n}/${d}`);

/** A polynomial with rational coefficients (constant first) as text. */
export function formatUnivariate(coeffs, v = "x") {
  const parts = [];
  for (let i = coeffs.length - 1; i >= 0; i -= 1) {
    const c = toQ(coeffs[i]);
    if (c[0] === 0n) continue;
    const neg = c[0] < 0n;
    const one = (c[0] === 1n || c[0] === -1n) && c[1] === 1n;
    const body = i === 0 ? fmtQ(c) : `${one ? "" : c[1] === 1n ? fmtQ(c) : `(${fmtQ(c)})`}${v}${i > 1 ? sup(i) : ""}`;
    parts.push((parts.length ? (neg ? " − " : " + ") : neg ? "−" : "") + body);
  }
  return parts.join("") || "0";
}

// ---- elements

export class TowerElement {
  constructor(field, v) {
    this.field = field;
    /** the same element in the underlying number field Q(θ) */
    this.v = v;
  }
  add(o) { return new TowerElement(this.field, this.v.add(this.field.of(o).v)); }
  sub(o) { return new TowerElement(this.field, this.v.sub(this.field.of(o).v)); }
  mul(o) { return new TowerElement(this.field, this.v.mul(this.field.of(o).v)); }
  neg() { return new TowerElement(this.field, this.v.neg()); }
  inverse() {
    if (this.v.isZero()) throw new RangeError("0 has no inverse");
    return new TowerElement(this.field, this.v.inv());
  }
  div(o) { return this.mul(this.field.of(o).inverse()); }
  pow(k) {
    if (!Number.isInteger(k)) throw new RangeError("the power must be a whole number");
    return new TowerElement(this.field, this.v.pow(k));
  }
  equals(o) { return this.v.equals(this.field.of(o).v); }
  isZero() { return this.v.isZero(); }
  /** Whether the element is a rational number. */
  get isRational() { return this.v.c.length <= 1; }
  /** The rational value, or null. */
  toRational() { return this.v.toRational(); }

  /** Coordinates in the tower basis (see `basis`), as Core Rationals. */
  coordinates() {
    return this.field._toTower(this.v).map(pub);
  }

  /** The element as a polynomial in the generators. */
  toPoly() {
    const K = this.field;
    const coords = K._toTower(this.v);
    return K.ring.fromTerms(coords.map((c, i) => [pub(c), K._exponents[i]]).filter(([c]) => c.numerator !== 0n));
  }

  /** The minimal polynomial over Q, constant first, as Core Rationals. */
  minimalPolynomial() {
    return this.field._minimal(this.v).map(pub);
  }
  /** The same as text in the variable x. */
  minimalPolynomialText(v = "x") { return formatUnivariate(this.field._minimal(this.v), v); }

  /** The degree of the element over Q: the degree of Q(this). */
  get degree() { return this.field._minimal(this.v).length - 1; }

  /**
   * The minimal polynomial over the field of the first k generators, constant
   * first, with coefficients in that field.
   */
  minimalPolynomialOver(k) {
    return this.field._minimalOver(this.v, k);
  }

  /** The norm and trace down to Q. */
  norm() { return this.v.norm(); }
  trace() { return this.v.trace(); }

  /** The exact value as an algebraic number, under the field's embedding into the complex numbers. */
  get value() { return this.v.toAlgebraic(); }
  /** A floating-point value: a number when real, { re, im } otherwise. */
  approx() {
    const z = this.value;
    if (z.rational) return Number(z.rational[0]) / Number(z.rational[1]);
    if (typeof z.toNumbers === "function") { const w = z.toNumbers(); return Array.isArray(w) ? { re: w[0], im: w[1] } : { re: w.re, im: w.im }; }
    return z.toNumber();
  }

  toString() { return this.toPoly().toString(); }
  toText() { return this.toPoly().toText(); }
  toTeX() { return this.toPoly().toTeX(); }
}

// ---- the tower

export class Extension {
  #L;
  #images = [];
  #toTower = null;
  #roots = null;
  #autos = null;

  constructor(steps) {
    if (!steps.length) throw new RangeError("an extension needs at least one generator");
    this.names = steps.map(([name]) => String(name));
    if (new Set(this.names).size !== this.names.length) throw new RangeError("the generators need different names");
    /** the polynomial ring in the generators, in which elements are written */
    this.ring = ring(this.names, { order: "lex" });
    const R = this.ring;
    let L = new NumberField([0n, 1n], { name: "θ", check: false });
    const degrees = [];
    const polys = [];
    steps.forEach(([name, text], k) => {
      const f = typeof text === "string" ? R.poly(text) : R.from(text);
      const later = f.variables().filter((v) => this.names.indexOf(v) > k);
      if (later.length) throw new RangeError(`the polynomial for ${name} uses ${later.join(", ")}, which comes later in the tower`);
      const coeffs = f.coefficients(name).map((c) => this.#evaluate(c, L));
      if (coeffs.length < 2) throw new RangeError(`the polynomial for ${name} must contain ${name}`);
      const { factors } = factorOver(L, coeffs);
      if (factors.length !== 1 || factors[0].multiplicity !== 1) {
        const over = k === 0 ? "Q" : `Q(${this.names.slice(0, k).join(", ")})`;
        const shape = factors.map((t) => `${t.factor.length - 1}${t.multiplicity > 1 ? `^${t.multiplicity}` : ""}`).join(", ");
        throw new RangeError(`${f} is not irreducible over ${over}: it has factors of degrees ${shape}, so ${name} would not generate a field of degree ${coeffs.length - 1} over it`);
      }
      const K = L;
      const ext = extend(K, coeffs, { name: "θ" });
      this.#images = this.#images.map((x) => ext.embed(x));
      this.#images.push(ext.root);
      if (ext.field !== K) chooseEmbedding(ext, K, k === 0);
      L = ext.field;
      degrees.push(coeffs.length - 1);
      polys.push(f);
    });
    this.#L = L;
    /** the degree of each step, and of the whole field over Q */
    this.degrees = degrees;
    this.degree = L.n;
    this.polynomials = polys;
    // the tower basis: products of powers of the generators, the first generator varying slowest
    const exps = [[]];
    for (const d of degrees) {
      const next = [];
      for (const e of exps) for (let i = 0; i < d; i += 1) next.push([...e, i]);
      exps.splice(0, exps.length, ...next);
    }
    this._exponents = exps;
    this.zero = new TowerElement(this, L.zero);
    this.one = new TowerElement(this, L.one);
    this.generators = this.#images.map((x) => new TowerElement(this, x));
    for (const [i, name] of this.names.entries()) if (!(name in this)) this[name] = this.generators[i];
  }

  #evaluate(poly, L = this.#L) {
    let out = L.zero;
    for (const { coefficient, exponents } of poly.terms()) {
      let t = L.of([coefficient.numerator, coefficient.denominator]);
      exponents.forEach((e, i) => { if (e) t = t.mul(this.#images[i].pow(e)); });
      out = out.add(t);
    }
    return out;
  }

  #vector(v) {
    return Array.from({ length: this.degree }, (_, i) => v.c[i] ?? Z);
  }

  /** Coordinates of a number field element in the tower basis. */
  _toTower(v) {
    if (!this.#toTower) {
      const cols = this._exponents.map((e) => {
        let t = this.#L.one;
        e.forEach((k, i) => { if (k) t = t.mul(this.#images[i].pow(k)); });
        return this.#vector(t);
      });
      this.#toTower = cols;
    }
    const x = solveColumns(this.#toTower, this.#vector(v));
    if (!x) throw new Error("tower: the basis does not span the field");
    return x;
  }

  _minimal(v) {
    const L = this.#L;
    const powers = [this.#vector(L.one)];
    let p = L.one;
    for (;;) {
      p = p.mul(v);
      const rel = solveColumns(powers, this.#vector(p));
      if (rel) return [...rel.map(([n, d]) => [-n, d]), ONE];
      powers.push(this.#vector(p));
    }
  }

  _minimalOver(v, k) {
    if (!Number.isInteger(k) || k < 0 || k > this.names.length) throw new RangeError(`the tower has ${this.names.length} generators`);
    const L = this.#L;
    const sub = this._exponents.filter((e) => e.slice(k).every((x) => x === 0));
    const base = sub.map((e) => { let t = L.one; e.forEach((x, i) => { if (x) t = t.mul(this.#images[i].pow(x)); }); return t; });
    const cols = [];
    let p = L.one;
    for (let d = 0; ; d += 1) {
      const rel = cols.length ? solveColumns(cols, this.#vector(p)) : null;
      if (rel) {
        const out = [];
        for (let j = 0; j < d; j += 1) {
          let c = L.zero;
          base.forEach((b, i) => { c = c.sub(b.mul(L.of(rel[j * base.length + i]))); });
          out.push(new TowerElement(this, c));
        }
        out.push(this.one);
        return out;
      }
      for (const b of base) cols.push(this.#vector(b.mul(p)));
      p = p.mul(v);
    }
  }

  /** An element from a number, a Rational, text in the generators such as "a*b + 1/(a + 1)", or a polynomial. */
  of(x) {
    if (x instanceof TowerElement) {
      if (x.field !== this) throw new TypeError("that element belongs to another field");
      return x;
    }
    if (typeof x === "string" && !/^\s*-?\d+\s*(\/\s*\d+\s*)?$/.test(x)) {
      if (x.includes("/")) {
        const q = fractions(this.ring).parse(x);
        return new TowerElement(this, this.#evaluate(q.num).div(this.#evaluate(q.den)));
      }
      return new TowerElement(this, this.#evaluate(this.ring.poly(x)));
    }
    if (x && typeof x === "object" && x.ring === this.ring) return new TowerElement(this, this.#evaluate(x));
    return new TowerElement(this, this.#L.of(toQ(x)));
  }

  /** The tower basis: the products of powers of the generators, as elements. */
  get basis() {
    return this._exponents.map((e) => new TowerElement(this, this.#evaluate(this.ring.monomial(e))));
  }

  /** The degree over the field of the first k generators. */
  degreeOver(k) {
    return this.degrees.slice(k).reduce((a, b) => a * b, 1);
  }

  /** The tower law written out: "[Q(a, b) : Q] = [Q(a, b) : Q(a)]·[Q(a) : Q] = 3·2 = 6". */
  towerLaw() {
    const name = (k) => (k === 0 ? "Q" : `Q(${this.names.slice(0, k).join(", ")})`);
    const r = this.names.length;
    const steps = [];
    for (let k = r; k >= 1; k -= 1) steps.push(`[${name(k)} : ${name(k - 1)}]`);
    if (r === 1) return `[${name(1)} : Q] = ${this.degree}`;
    return `[${name(r)} : Q] = ${steps.join("·")} = ${[...this.degrees].reverse().join("·")} = ${this.degree}`;
  }

  /**
   * A single generator of the whole field: { element, polynomial (its minimal
   * polynomial, constant first), text, expressions } where `expressions`
   * writes each tower generator as a polynomial in the primitive element θ.
   */
  primitiveElement() {
    const L = this.#L;
    const expressions = {};
    this.names.forEach((nm, i) => { expressions[nm] = formatUnivariate(this.#vector(this.#images[i]), "θ"); });
    return { element: new TowerElement(this, L.gen), polynomial: L.m.map(pub), text: formatUnivariate(L.m, "x"), expressions };
  }

  /** The field Q(x) generated by one element: { degree, polynomial, text, isWhole }. */
  subfield(x) {
    const m = this._minimal(this.of(x).v);
    return { degree: m.length - 1, polynomial: m.map(pub), text: formatUnivariate(m, "x"), isWhole: m.length - 1 === this.degree };
  }

  #polyOver(coeffs) {
    const list = typeof coeffs === "string" ? this.#readUnivariate(coeffs) : coeffs;
    return list.map((c) => this.of(c).v);
  }

  #readUnivariate(text) {
    // a polynomial in one extra variable x (or the first letter that is not a generator)
    const letters = [...new Set(String(text).match(/[A-Za-z_Ͱ-Ͽ]\w*/g) ?? [])].filter((s) => !this.names.includes(s));
    if (letters.length !== 1) throw new RangeError("give the polynomial in one variable besides the generators, such as \"x^2 - 6\"");
    const R = ring([letters[0], ...this.names], { order: "lex" });
    const inner = (p) => this.ring.fromTerms(p.terms().map((t) => [t.coefficient, t.exponents.slice(1)]));
    return R.poly(text).coefficients(letters[0]).map((c) => this.of(inner(c)));
  }

  /** The roots of a polynomial (text in x, or coefficients constant first) that lie in the field. */
  roots(coeffs) {
    return rootsIn(this.#L, this.#polyOver(coeffs)).map((r) => ({ value: new TowerElement(this, r.value), multiplicity: r.multiplicity }));
  }

  /** The irreducible factors of a polynomial over this field: [{ factor: [elements, constant first], multiplicity, degree, text }]. */
  factor(coeffs, { variable = "x" } = {}) {
    return factorOver(this.#L, this.#polyOver(coeffs)).factors.map((f) => {
      const els = f.factor.map((c) => new TowerElement(this, c));
      return { factor: els, multiplicity: f.multiplicity, degree: els.length - 1, text: formatOver(els, variable) };
    });
  }

  #conjugateRoots() {
    if (!this.#roots) {
      const L = this.#L;
      const found = rootsIn(L, L.m.map((c) => L.of(c))).map((r) => r.value);
      // the identity (θ itself) first
      found.sort((a, b) => (a.equals(L.gen) ? -1 : b.equals(L.gen) ? 1 : 0));
      this.#roots = found;
    }
    return this.#roots;
  }

  #applyRoot(rho, v) {
    const L = this.#L;
    let acc = L.zero;
    for (let i = v.c.length - 1; i >= 0; i -= 1) acc = acc.mul(rho).add(L.of(v.c[i]));
    return acc;
  }

  /**
   * The automorphisms of the field: each { images, apply(x), toString() },
   * the identity first. There are as many as the degree exactly when the
   * field is a Galois extension of Q.
   */
  automorphisms() {
    if (!this.#autos) {
      this.#autos = this.#conjugateRoots().map((rho, index) => {
        const images = {};
        this.names.forEach((nm, i) => { images[nm] = new TowerElement(this, this.#applyRoot(rho, this.#images[i])); });
        const text = this.names.map((nm) => `${nm} ↦ ${images[nm]}`).join(", ");
        return Object.freeze({ index, images, apply: (x) => new TowerElement(this, this.#applyRoot(rho, this.of(x).v)), toString: () => text });
      });
    }
    return this.#autos;
  }

  /** Whether the field is a Galois extension of Q (normal: every conjugate of every element lies in it). */
  get isGalois() { return this.#conjugateRoots().length === this.degree; }

  /**
   * The group of automorphisms as a `Group`, with elements the automorphisms
   * (multiplied as maps: first the right one). For a Galois extension this is
   * the Galois group.
   */
  automorphismGroup({ name = null } = {}) {
    const autos = this.automorphisms();
    const rhos = this.#conjugateRoots();
    const mul = (s, t) => {
      const image = this.#applyRoot(rhos[s.index], rhos[t.index]);
      return autos[rhos.findIndex((r) => r.equals(image))];
    };
    return Group.fromElements(autos, mul, { key: (s) => String(s.index), format: (s) => (s.index === 0 ? "id" : String(s)), identity: autos[0], name: name ?? (this.isGalois ? `Gal(${this}/Q)` : `Aut(${this})`) });
  }
  galoisGroup(options) {
    if (!this.isGalois) throw new RangeError(`${this} is not a Galois extension of Q: it has ${this.automorphisms().length} automorphisms and degree ${this.degree}`);
    return this.automorphismGroup(options);
  }

  /**
   * The elements fixed by some automorphisms: { degree, basis, generator,
   * polynomial, text }. `generator` is one element generating the fixed
   * field and `polynomial` its minimal polynomial.
   */
  fixedField(automorphisms) {
    const list = (Array.isArray(automorphisms) ? automorphisms : [automorphisms]).map((s) => (typeof s === "number" ? this.automorphisms()[s] : s));
    const n = this.degree;
    const L = this.#L;
    const rhos = this.#conjugateRoots();
    const unit = (i) => L.of(Array.from({ length: n }, (_, j) => (j === i ? ONE : Z)));
    // rows of (σ − 1) for each σ, acting on coordinates in the powers of θ
    const rows = [];
    for (const s of list) {
      const cols = Array.from({ length: n }, (_, i) => this.#vector(this.#applyRoot(rhos[s.index], unit(i)).sub(unit(i))));
      for (let r = 0; r < n; r += 1) rows.push(cols.map((c) => c[r]));
    }
    const kernel = rows.length ? nullSpace(rows, n) : Array.from({ length: n }, (_, i) => Array.from({ length: n }, (_, j) => (j === i ? ONE : Z)));
    const basis = kernel.map((v) => new TowerElement(this, L.of(v)));
    const dim = basis.length;
    // a generator: first the sums over the orbit of each basis monomial, then small combinations of the basis
    let generator = null;
    let minimal = null;
    const accept = (el) => {
      if (generator) return;
      const m = this._minimal(el);
      if (m.length - 1 === dim) { generator = new TowerElement(this, el); minimal = m; }
    };
    if (dim === 1) accept(L.one);
    if (dim === n) accept(L.gen);
    const closure = this.#closure(list.map((s) => s.index));
    const order = this._exponents.map((e, i) => i).sort((p, q) => this._exponents[p].reduce((a, b) => a + b, 0) - this._exponents[q].reduce((a, b) => a + b, 0));
    for (const i of order) {
      if (generator) break;
      const mono = this.#evaluate(this.ring.monomial(this._exponents[i]));
      const seen = [];
      for (const k of closure) { const img = this.#applyRoot(rhos[k], mono); if (!seen.some((x) => x.equals(img))) seen.push(img); }
      accept(seen.reduce((a, b) => a.add(b), L.zero));
      if (seen.length > 1) accept(seen.reduce((a, b) => a.mul(b), L.one));
    }
    for (let bound = 1; bound <= 3 && !generator; bound += 1) {
      const width = Math.min(dim, 5);
      const total = (2 * bound + 1) ** width;
      for (let t = 0; t < total && !generator; t += 1) {
        let x = t;
        let el = L.zero;
        for (let i = 0; i < width; i += 1) { const c = (x % (2 * bound + 1)) - bound; x = Math.floor(x / (2 * bound + 1)); if (c) el = el.add(basis[i].v.mul(L.of(BigInt(c)))); }
        accept(el);
      }
    }
    return { degree: dim, basis, generator, polynomial: minimal ? minimal.map(pub) : null, text: minimal ? formatUnivariate(minimal, "x") : null };
  }

  /** The indices of the automorphisms generated by the given ones. */
  #closure(indices) {
    const rhos = this.#conjugateRoots();
    const seen = new Set([0, ...indices]);
    const queue = [...seen];
    for (let h = 0; h < queue.length; h += 1) for (const g of indices) {
      const image = this.#applyRoot(rhos[g], rhos[queue[h]]);
      const k = rhos.findIndex((r) => r.equals(image));
      if (!seen.has(k)) { seen.add(k); queue.push(k); }
    }
    return queue;
  }

  /** The field as a NumberField of @ratmath/algebraic, with the generators' images in it. */
  toNumberField() {
    return { field: this.#L, generators: [...this.#images] };
  }

  toString() { return `Q(${this.names.join(", ")})`; }
  describe() {
    // each polynomial printed with its own generator leading
    const shown = this.polynomials.map((f, k) => f.into(ring([this.names[k], ...this.names.filter((_, i) => i !== k)], { order: "lex" })));
    return `${this}: degree ${this.degree} over Q, with ${shown.map((f) => `${f} = 0`).join(", ")}`;
  }
}

function formatOver(els, v) {
  const parts = [];
  for (let i = els.length - 1; i >= 0; i -= 1) {
    if (els[i].isZero()) continue;
    const r = els[i].toRational();
    let c;
    let neg = false;
    if (r) { neg = r.numerator < 0n; const a = neg ? r.negate() : r; c = i > 0 && a.numerator === 1n && a.denominator === 1n ? "" : a.denominator === 1n ? String(a) : `(${a})`; }
    else {
      const t = String(els[i]);
      if (/ [+−] /.test(t)) c = `(${t})`;
      else { neg = t.startsWith("−"); c = neg ? t.slice(1) : t; if (i > 0 && /^\d/.test(c) === false) c = `${c}·`; else if (i > 0) c = `${c}·`; }
    }
    const body = i === 0 ? c || "1" : `${c}${v}${i > 1 ? sup(i) : ""}`;
    parts.push((parts.length ? (neg ? " − " : " + ") : neg ? "−" : "") + body);
  }
  return parts.join("") || "0";
}

/**
 * Chooses which complex number the new primitive element stands for: one
 * that keeps the earlier generators where they were and makes the new
 * generator its largest real root when it has a real root, else a root in
 * the upper half plane.
 */
function chooseEmbedding(ext, K, first) {
  const L = ext.field;
  const alpha = first ? null : ext.embed(K.gen);
  const beta = ext.root;
  const evalAt = (el, [re, im]) => {
    let a = 0;
    let b = 0;
    for (let i = el.c.length - 1; i >= 0; i -= 1) {
      const c = Number(el.c[i][0]) / Number(el.c[i][1]);
      [a, b] = [a * re - b * im + c, a * im + b * re];
    }
    return [a, b];
  };
  const numbers = (z) => {
    if (z.rational) return [Number(z.rational[0]) / Number(z.rational[1]), 0];
    if (typeof z.toNumbers !== "function") return [z.toNumber(), 0];
    const w = z.toNumbers();
    return Array.isArray(w) ? w : [w.re, w.im];
  };
  const previous = L.rootChoice;
  L.rootChoice = () => {
    const all = algebraicRoots(L.m, { complex: true }).map((r) => r.value);
    const scale = (z) => 1 + Math.hypot(z[0], z[1]);
    let candidates = all.map((r) => ({ r, z: numbers(r) }));
    if (alpha) {
      const target = numbers(K.root);
      const close = candidates.filter(({ z }) => { const w = evalAt(alpha, z); return Math.hypot(w[0] - target[0], w[1] - target[1]) < 1e-7 * scale(target); });
      if (close.length !== candidates.length / K.n) return typeof previous === "function" ? previous() : all[0];
      candidates = close;
    }
    const scored = candidates.map(({ r, z }) => ({ r, w: evalAt(beta, z) }));
    const real = scored.filter(({ w }) => Math.abs(w[1]) < 1e-9 * scale(w)).sort((p, q) => q.w[0] - p.w[0]);
    if (real.length) return real[0].r;
    const upper = scored.filter(({ w }) => w[1] > 0).sort((p, q) => q.w[1] - p.w[1]);
    return (upper[0] ?? scored[0]).r;
  };
}

/**
 * A tower of extensions of Q: extension({ a: "a^2 - 2", b: "b^3 - a" }), or
 * a list of [name, polynomial] pairs. Each polynomial is in its own
 * generator, with coefficients that may use the generators before it.
 */
export function extension(steps) {
  const list = Array.isArray(steps) ? steps : Object.entries(steps);
  return new Extension(list.map((s) => (Array.isArray(s) ? s : [s.name, s.polynomial])));
}
