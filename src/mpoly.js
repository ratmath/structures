// Polynomials in several variables with exact coefficients.
//
//   const R = ring("x, y, z");                 // Q[x, y, z], graded reverse lex
//   const f = R.poly("x^2 y - 3/2 x z + 1");
//   f.mul(f).degree();                         // 6
//
// A polynomial is a sparse sum of terms: a map from exponent vectors to
// non-zero coefficients in the ring's field. Polynomials are immutable.

import { Rational } from "@ratmath/core";
import { QQ, normalize, fromInteger } from "./field.js";
import { orderOf, keyOf, mulExp, dividesExp, divExp, totalDegree, superscript } from "./monomial.js";
import { parsePolynomial } from "./parse.js";

const NAME = /^[A-Za-z][A-Za-z0-9_]*$/;

export class PolyRing {
  /**
   * @param vars    variable names: an array, or text such as "x, y, z" or "x y z"
   * @param options { field: a field object (default the rationals), order: "grevlex" | "grlex" | "lex" | function }
   */
  constructor(vars, { field, order = "grevlex" } = {}) {
    const names = typeof vars === "string" ? vars.split(/[\s,]+/).filter(Boolean) : [...vars];
    if (!names.length) throw new RangeError("a polynomial ring needs at least one variable");
    for (const v of names) {
      if (typeof v !== "string" || !NAME.test(v)) throw new RangeError(`"${v}" is not a variable name (a letter, then letters, digits or _)`);
    }
    if (new Set(names).size !== names.length) throw new RangeError("variable names must be distinct");
    this.vars = Object.freeze(names);
    this.n = names.length;
    this.K = normalize(field);
    this.field = field ?? QQ;
    this.order = typeof order === "function" ? "custom" : order;
    this.cmp = orderOf(order);
    this._index = new Map(names.map((v, i) => [v, i]));
    this._zeroExp = Object.freeze(new Array(this.n).fill(0));
    this.zero = new MPoly(this, new Map());
    this.one = this.constant(this.K.one);
    this.gens = Object.freeze(names.map((_, i) => this._gen(i)));
    Object.freeze(this);
  }

  _gen(i) {
    const e = new Array(this.n).fill(0);
    e[i] = 1;
    return new MPoly(this, new Map([[keyOf(e), { e, c: this.K.one }]]));
  }

  /** The position of a variable given by name or position. */
  indexOf(v) {
    if (typeof v === "number") {
      if (!Number.isInteger(v) || v < 0 || v >= this.n) throw new RangeError(`there is no variable number ${v}`);
      return v;
    }
    if (v instanceof MPoly) {
      const i = v._asVariable();
      if (i >= 0 && v.ring === this) return i;
      throw new TypeError("expected a variable");
    }
    const i = this._index.get(v);
    if (i === undefined) throw new RangeError(`"${v}" is not a variable of ${this}`);
    return i;
  }

  /** The variable with this name, as a polynomial. */
  gen(v) {
    return this.gens[this.indexOf(v)];
  }

  /** A constant polynomial from a field element (already in the field's own representation). */
  constant(c) {
    return this.K.sign(c) === 0 ? this.zero : new MPoly(this, new Map([[keyOf(this._zeroExp), { e: this._zeroExp, c }]]));
  }

  /** c·x^exps, with exps an array of exponents or an object such as { x: 2, y: 1 }. */
  monomial(exps, c = 1) {
    const e = this._exps(exps);
    const coeff = this._coeff(c);
    return this.K.sign(coeff) === 0 ? this.zero : new MPoly(this, new Map([[keyOf(e), { e, c: coeff }]]));
  }

  _exps(exps) {
    let e;
    if (Array.isArray(exps)) {
      if (exps.length !== this.n) throw new RangeError(`expected ${this.n} exponents`);
      e = [...exps];
    } else {
      e = new Array(this.n).fill(0);
      for (const [v, k] of Object.entries(exps)) e[this.indexOf(v)] = k;
    }
    for (const k of e) if (!Number.isInteger(k) || k < 0) throw new RangeError("exponents must be non-negative integers");
    return e;
  }

  // a coefficient from caller input: numbers, bigints, text, Core Rationals, or an element of the field
  _coeff(x) {
    const K = this.K;
    if (K === QQ) return QQ.from(x);
    if (typeof x === "number" || typeof x === "bigint" || typeof x === "string" || x instanceof Rational) return K.from(x);
    return x;
  }

  /**
   * Anything reasonable as a polynomial of this ring: a polynomial of this
   * ring, text ("x^2 - y"), or a coefficient.
   */
  from(x) {
    if (x instanceof MPoly) {
      if (x.ring === this) return x;
      return x.into(this);
    }
    if (typeof x === "string") return this.poly(x);
    return this.constant(this._coeff(x));
  }

  /** Read a polynomial from text: + - * ^, parentheses, exact numbers, division by constants. */
  poly(text) {
    return parsePolynomial(this, String(text));
  }

  /** Build from a list of [coefficient, exponents] pairs. */
  fromTerms(terms) {
    const m = new Map();
    for (const [c, exps] of terms) addTerm(this.K, m, this._exps(exps), this._coeff(c));
    return new MPoly(this, m);
  }

  /**
   * A polynomial in one variable from its coefficients, constant first (the
   * form @ratmath/algebraic uses).
   */
  fromUnivariate(coeffs, v = 0) {
    const i = this.indexOf(v);
    const m = new Map();
    coeffs.forEach((c, k) => {
      const e = new Array(this.n).fill(0);
      e[i] = k;
      addTerm(this.K, m, e, this._coeff(c));
    });
    return new MPoly(this, m);
  }

  /** The same variables and field with another monomial order. */
  withOrder(order) {
    return new PolyRing(this.vars, { field: this.field, order });
  }

  equals(other) {
    return other === this || (other instanceof PolyRing && other.field === this.field && other.order === this.order && other.cmp === this.cmp && other.n === this.n && other.vars.every((v, i) => v === this.vars[i]));
  }

  toString() {
    return `${this.K.name}[${this.vars.join(", ")}]`;
  }
}

/** A polynomial ring: `ring("x, y")`, `ring(["x", "y"], { field: primeField(7), order: "lex" })`. */
export const ring = (vars, options) => new PolyRing(vars, options);

function addTerm(K, m, e, c) {
  if (K.sign(c) === 0) return;
  const k = keyOf(e);
  const old = m.get(k);
  if (!old) {
    m.set(k, { e, c });
    return;
  }
  const s = K.add(old.c, c);
  if (K.sign(s) === 0) m.delete(k);
  else m.set(k, { e, c: s });
}

export class MPoly {
  constructor(ringOf, terms) {
    this.ring = ringOf;
    this._m = terms;
    this._sorted = null;
  }

  // ---- inspection ---------------------------------------------------------

  /** Terms as { coefficient, exponents }, largest first in the ring's order. */
  _terms() {
    if (!this._sorted) {
      const cmp = this.ring.cmp;
      this._sorted = [...this._m.values()].sort((a, b) => cmp(b.e, a.e));
    }
    return this._sorted;
  }

  /** Terms, largest first: [{ coefficient, exponents, monomial }]. */
  terms() {
    const { ring: R } = this;
    return this._terms().map(({ e, c }) => ({ coefficient: R.K.toPublic(c), exponents: [...e], monomial: new MPoly(R, new Map([[keyOf(e), { e, c: R.K.one }]])) }));
  }

  get length() {
    return this._m.size;
  }

  isZero() {
    return this._m.size === 0;
  }

  isConstant() {
    return this._m.size === 0 || (this._m.size === 1 && this._m.has(keyOf(this.ring._zeroExp)));
  }

  isMonomial() {
    return this._m.size === 1;
  }

  // the index of the variable this polynomial is, or −1
  _asVariable() {
    if (this._m.size !== 1) return -1;
    const [{ e, c }] = this._m.values();
    if (!this.ring.K.equals(c, this.ring.K.one) || totalDegree(e) !== 1) return -1;
    return e.indexOf(1);
  }

  /** The constant term, as a field element. */
  constantTerm() {
    const t = this._m.get(keyOf(this.ring._zeroExp));
    return this.ring.K.toPublic(t ? t.c : this.ring.K.zero);
  }

  /** The value of a constant polynomial; throws if it is not constant. */
  toConstant() {
    if (!this.isConstant()) throw new RangeError(`${this} is not a constant`);
    return this.constantTerm();
  }

  /**
   * Total degree, or the degree in one variable. The zero polynomial has
   * degree −Infinity.
   */
  degree(v) {
    if (this.isZero()) return -Infinity;
    let d = 0;
    if (v === undefined) {
      for (const { e } of this._m.values()) d = Math.max(d, totalDegree(e));
    } else {
      const i = this.ring.indexOf(v);
      for (const { e } of this._m.values()) d = Math.max(d, e[i]);
    }
    return d;
  }

  /**
   * The coefficients, constant first, of a polynomial that involves at most
   * one variable (the form @ratmath/algebraic uses). Throws otherwise.
   */
  toUnivariate() {
    const used = this.variables();
    if (used.length > 1) throw new RangeError(`${this} involves more than one variable`);
    const K = this.ring.K;
    if (!used.length) return this.isZero() ? [] : [K.toPublic(this._lead().c)];
    return this.coefficients(used[0]).map((c) => K.toPublic(c.isZero() ? K.zero : c._lead().c));
  }

  /** The degree in each variable, as an object keyed by variable name. */
  degrees() {
    return Object.fromEntries(this.ring.vars.map((v) => [v, this.degree(v)]));
  }

  /** The names of the variables that actually occur. */
  variables() {
    const used = new Array(this.ring.n).fill(false);
    for (const { e } of this._m.values()) e.forEach((k, i) => k > 0 && (used[i] = true));
    return this.ring.vars.filter((_, i) => used[i]);
  }

  isHomogeneous() {
    let d = null;
    for (const { e } of this._m.values()) {
      const t = totalDegree(e);
      if (d === null) d = t;
      else if (d !== t) return false;
    }
    return true;
  }

  _lead() {
    if (this.isZero()) throw new RangeError("the zero polynomial has no leading term");
    return this._terms()[0];
  }

  /** The largest term in the ring's order, as a polynomial. */
  leadingTerm() {
    const { e, c } = this._lead();
    return new MPoly(this.ring, new Map([[keyOf(e), { e, c }]]));
  }

  leadingCoefficient() {
    return this.ring.K.toPublic(this._lead().c);
  }

  leadingMonomial() {
    const { e } = this._lead();
    return new MPoly(this.ring, new Map([[keyOf(e), { e, c: this.ring.K.one }]]));
  }

  leadingExponents() {
    return [...this._lead().e];
  }

  /** The coefficient of a monomial given as exponents ([2, 1] or { x: 2, y: 1 }). */
  coefficient(exps) {
    const t = this._m.get(keyOf(this.ring._exps(exps)));
    return this.ring.K.toPublic(t ? t.c : this.ring.K.zero);
  }

  /**
   * As a polynomial in one variable: the list of coefficients, constant first,
   * each a polynomial in the other variables.
   */
  coefficients(v) {
    const i = this.ring.indexOf(v);
    const out = [];
    for (const { e, c } of this._m.values()) {
      const k = e[i];
      while (out.length <= k) out.push(new Map());
      const rest = [...e];
      rest[i] = 0;
      out[k].set(keyOf(rest), { e: rest, c });
    }
    return out.map((m) => new MPoly(this.ring, m));
  }

  // ---- arithmetic ---------------------------------------------------------

  _same(other) {
    const g = this.ring.from(other);
    return g;
  }

  add(other) {
    const g = this._same(other);
    if (g.isZero()) return this;
    if (this.isZero()) return g;
    const m = new Map(this._m);
    const K = this.ring.K;
    for (const { e, c } of g._m.values()) addTerm(K, m, e, c);
    return new MPoly(this.ring, m);
  }

  neg() {
    const K = this.ring.K;
    const m = new Map();
    for (const [k, { e, c }] of this._m) m.set(k, { e, c: K.neg(c) });
    return new MPoly(this.ring, m);
  }

  sub(other) {
    return this.add(this._same(other).neg());
  }

  mul(other) {
    const g = this._same(other);
    if (this.isZero() || g.isZero()) return this.ring.zero;
    const K = this.ring.K;
    const m = new Map();
    for (const a of this._m.values()) {
      for (const b of g._m.values()) addTerm(K, m, mulExp(a.e, b.e), K.mul(a.c, b.c));
    }
    return new MPoly(this.ring, m);
  }

  /** Multiply by a coefficient. */
  scale(c) {
    return this._scale(this.ring._coeff(c));
  }

  _scale(c) {
    const K = this.ring.K;
    if (K.sign(c) === 0) return this.ring.zero;
    const m = new Map();
    for (const [k, t] of this._m) {
      const p = K.mul(t.c, c);
      if (K.sign(p) !== 0) m.set(k, { e: t.e, c: p });
    }
    return new MPoly(this.ring, m);
  }

  // multiply by the term c·x^e
  _mulTerm(e, c) {
    const K = this.ring.K;
    const m = new Map();
    for (const t of this._m.values()) {
      const ee = mulExp(t.e, e);
      const p = K.mul(t.c, c);
      if (K.sign(p) !== 0) m.set(keyOf(ee), { e: ee, c: p });
    }
    return new MPoly(this.ring, m);
  }

  pow(n) {
    if (typeof n === "bigint") n = Number(n);
    if (!Number.isInteger(n) || n < 0) throw new RangeError("a polynomial can only be raised to a non-negative integer power");
    let result = this.ring.one;
    let base = this;
    while (n > 0) {
      if (n & 1) result = result.mul(base);
      n >>= 1;
      if (n) base = base.mul(base);
    }
    return result;
  }

  equals(other) {
    const g = this._same(other);
    if (g._m.size !== this._m.size) return false;
    const K = this.ring.K;
    for (const [k, { c }] of this._m) {
      const t = g._m.get(k);
      if (!t || !K.equals(t.c, c)) return false;
    }
    return true;
  }

  /** Divide every coefficient by the leading one, so the leading coefficient is 1. */
  monic() {
    if (this.isZero()) return this;
    return this._scale(this.ring.K.inv(this._lead().c));
  }

  // ---- division -----------------------------------------------------------

  /**
   * The division algorithm in several variables. Divides by one polynomial or
   * a list of them, in the ring's monomial order:
   *
   *   f = q₁·g₁ + … + qₖ·gₖ + r
   *
   * where no term of r is divisible by the leading term of any gᵢ. With one
   * divisor, r is zero exactly when it divides f. With several, the remainder
   * can depend on their order unless they form a Gröbner basis.
   */
  divmod(divisors) {
    const R = this.ring;
    const K = R.K;
    const single = !Array.isArray(divisors);
    const gs = (single ? [divisors] : divisors).map((g) => R.from(g));
    if (!gs.length) throw new RangeError("nothing to divide by");
    const leads = gs.map((g) => {
      if (g.isZero()) throw new RangeError("division by the zero polynomial");
      return g._lead();
    });
    const qs = gs.map(() => new Map());
    const r = new Map();
    // working copy of f as a map; repeatedly remove its leading term
    const p = new Map(this._m);
    const cmp = R.cmp;
    while (p.size) {
      let top = null;
      for (const t of p.values()) if (!top || cmp(t.e, top.e) > 0) top = t;
      let divided = false;
      for (let i = 0; i < gs.length; i += 1) {
        const lt = leads[i];
        if (!dividesExp(lt.e, top.e)) continue;
        const e = divExp(top.e, lt.e);
        const c = K.div(top.c, lt.c);
        addTerm(K, qs[i], e, c);
        const nc = K.neg(c);
        for (const t of gs[i]._m.values()) addTerm(K, p, mulExp(t.e, e), K.mul(t.c, nc));
        divided = true;
        break;
      }
      if (!divided) {
        addTerm(K, r, top.e, top.c);
        p.delete(keyOf(top.e));
      }
    }
    const quotients = qs.map((m) => new MPoly(R, m));
    const remainder = new MPoly(R, r);
    return single ? { quotient: quotients[0], remainder } : { quotients, remainder };
  }

  /** The remainder on division by one polynomial or a list (see divmod). */
  rem(divisors) {
    return this.divmod(divisors).remainder;
  }

  /** Whether g divides this polynomial exactly. */
  isDivisibleBy(g) {
    return this.divmod(this.ring.from(g)).remainder.isZero();
  }

  /** The exact quotient; throws if the division leaves a remainder. */
  div(g) {
    g = this.ring.from(g);
    if (g.isConstant()) {
      if (g.isZero()) throw new RangeError("division by the zero polynomial");
      return this._scale(this.ring.K.inv(g._lead().c));
    }
    const { quotient, remainder } = this.divmod(g);
    if (!remainder.isZero()) throw new RangeError(`${g} does not divide ${this} exactly`);
    return quotient;
  }

  // ---- calculus and substitution -----------------------------------------

  /** The partial derivative with respect to a variable (n times). */
  diff(v, times = 1) {
    const i = this.ring.indexOf(v);
    const K = this.ring.K;
    let f = this;
    for (let t = 0; t < times; t += 1) {
      const m = new Map();
      for (const { e, c } of f._m.values()) {
        if (e[i] === 0) continue;
        const ee = [...e];
        ee[i] -= 1;
        const cc = K.mul(c, fromInteger(K, e[i]));
        if (K.sign(cc) !== 0) m.set(keyOf(ee), { e: ee, c: cc });
      }
      f = new MPoly(this.ring, m);
    }
    return f;
  }

  /** The gradient: the list of partial derivatives, in the ring's variable order. */
  gradient() {
    return this.ring.vars.map((v) => this.diff(v));
  }

  /**
   * Substitute for some of the variables. Values may be coefficients or
   * polynomials of this ring; variables not mentioned stay. The result is a
   * polynomial of the same ring: `f.subs({ x: 2 })`, `f.subs({ x: "y + 1" })`.
   */
  subs(values) {
    const R = this.ring;
    const repl = new Array(R.n).fill(null);
    for (const [v, x] of Object.entries(values)) repl[R.indexOf(v)] = R.from(x);
    if (repl.every((x) => x === null)) return this;
    // powers of each replacement, computed as needed
    const pows = repl.map((x) => (x ? [R.one, x] : null));
    const power = (i, k) => {
      const list = pows[i];
      while (list.length <= k) list.push(list[list.length - 1].mul(repl[i]));
      return list[k];
    };
    let out = R.zero;
    for (const { e, c } of this._m.values()) {
      const rest = [...e];
      let term = null;
      for (let i = 0; i < R.n; i += 1) {
        if (!repl[i] || e[i] === 0) continue;
        rest[i] = 0;
        term = term ? term.mul(power(i, e[i])) : power(i, e[i]);
      }
      const base = new MPoly(R, new Map([[keyOf(rest), { e: rest, c }]]));
      out = out.add(term ? base.mul(term) : base);
    }
    return out;
  }

  /**
   * Replace every variable by a polynomial of another ring (a change of
   * variables, or evaluating at polynomials): `f.compose({ s: "x + y", p: "x y" }, R)`.
   * Every variable that occurs must be given. The target ring is `target`, or
   * the ring of the first polynomial among the values.
   */
  compose(values, target) {
    const R = this.ring;
    const list = Array.isArray(values) ? values : R.vars.map((v) => values[v]);
    if (Array.isArray(values) && values.length !== R.n) throw new RangeError(`expected ${R.n} polynomials`);
    const T = target ?? list.find((x) => x instanceof MPoly)?.ring;
    if (!T) throw new TypeError("compose needs a target ring: pass one, or give polynomials");
    const repl = list.map((x) => (x === undefined || x === null ? null : T.from(x)));
    const pows = repl.map((x) => (x ? [T.one, x] : null));
    let out = T.zero;
    for (const { e, c } of this._m.values()) {
      let term = T.constant(T.K === R.K ? c : T.K.from(R.K.toPublic(c)));
      for (let i = 0; i < R.n; i += 1) {
        if (e[i] === 0) continue;
        if (!repl[i]) throw new RangeError(`no polynomial was given for ${R.vars[i]}`);
        const l = pows[i];
        while (l.length <= e[i]) l.push(l[l.length - 1].mul(repl[i]));
        term = term.mul(l[e[i]]);
      }
      out = out.add(term);
    }
    return out;
  }

  /**
   * The value at a point: every variable that occurs must be given. Returns a
   * field element (a Core Rational over Q).
   */
  at(point) {
    const R = this.ring;
    const K = R.K;
    const vals = new Array(R.n).fill(null);
    if (Array.isArray(point)) {
      if (point.length !== R.n) throw new RangeError(`expected ${R.n} values`);
      point.forEach((x, i) => (vals[i] = R._coeff(x)));
    } else {
      for (const [v, x] of Object.entries(point)) vals[R.indexOf(v)] = R._coeff(x);
    }
    const pows = vals.map((x) => (x === null ? null : [K.one, x]));
    let sum = K.zero;
    for (const { e, c } of this._m.values()) {
      let term = c;
      for (let i = 0; i < R.n; i += 1) {
        if (e[i] === 0) continue;
        if (!pows[i]) throw new RangeError(`no value was given for ${R.vars[i]}`);
        const list = pows[i];
        while (list.length <= e[i]) list.push(K.mul(list[list.length - 1], vals[i]));
        term = K.mul(term, list[e[i]]);
      }
      sum = K.add(sum, term);
    }
    return K.toPublic(sum);
  }

  /** Make homogeneous with a new last variable (the ring gains one variable). */
  homogenize(name = "h") {
    const R = new PolyRing([...this.ring.vars, name], { field: this.ring.field, order: this.ring.order === "custom" ? this.ring.cmp : this.ring.order });
    const d = this.isZero() ? 0 : this.degree();
    const m = new Map();
    for (const { e, c } of this._m.values()) {
      const ee = [...e, d - totalDegree(e)];
      m.set(keyOf(ee), { e: ee, c });
    }
    return new MPoly(R, m);
  }

  /**
   * The same polynomial in another ring. The target must contain every
   * variable that occurs; its field must accept this ring's coefficients
   * (the same field, or from Q into any field).
   */
  into(target) {
    const R = this.ring;
    if (target === R) return this;
    const map = R.vars.map((v) => target._index.get(v));
    const sameField = target.field === R.field;
    const m = new Map();
    for (const { e, c } of this._m.values()) {
      const ee = new Array(target.n).fill(0);
      for (let i = 0; i < R.n; i += 1) {
        if (e[i] === 0) continue;
        if (map[i] === undefined) throw new RangeError(`${target} has no variable ${R.vars[i]}`);
        ee[map[i]] = e[i];
      }
      let cc = c;
      if (!sameField) {
        if (R.K !== QQ) throw new TypeError(`cannot move coefficients from ${R.K.name} to ${target.K.name}`);
        cc = target.K.from(QQ.toPublic(c));
      }
      addTerm(target.K, m, ee, cc);
    }
    return new MPoly(target, m);
  }

  // ---- printing -----------------------------------------------------------

  /** Unicode text, largest term first: "x²y − (3/2)xz + 1". */
  toString() {
    return format(this, "unicode");
  }

  /** Plain text that `ring.poly` reads back: "x^2*y - 3/2*x*z + 1". */
  toText() {
    return format(this, "text");
  }

  toTeX() {
    return format(this, "tex");
  }

  toJSON() {
    return { vars: [...this.ring.vars], terms: this._terms().map(({ e, c }) => [this.ring.K.format(c), [...e]]) };
  }
}

function monomialText(R, e, style) {
  const parts = [];
  for (let i = 0; i < R.n; i += 1) {
    if (e[i] === 0) continue;
    const v = R.vars[i];
    if (style === "tex") {
      const name = v.length > 1 ? (/^[A-Za-z]\d+$/.test(v) ? `${v[0]}_{${v.slice(1)}}` : `\\mathit{${v.replace(/_/g, "\\_")}}`) : v;
      parts.push(e[i] === 1 ? name : `${name}^{${e[i]}}`);
    } else if (style === "unicode") parts.push(e[i] === 1 ? v : `${v}${superscript(e[i])}`);
    else parts.push(e[i] === 1 ? v : `${v}^${e[i]}`);
  }
  if (style === "text") return parts.join("*");
  // juxtaposition is only unambiguous for one-letter names
  const plain = R.vars.every((v) => v.length === 1);
  return parts.join(style === "tex" ? (plain ? "" : "\\,") : plain ? "" : "·");
}

function format(f, style) {
  const R = f.ring;
  const K = R.K;
  if (f.isZero()) return "0";
  const minus = style === "unicode" ? "−" : "-";
  const out = [];
  f._terms().forEach(({ e, c }, idx) => {
    const mono = monomialText(R, e, style);
    // work from the coefficient's own text, so every field prints the same way
    let text = K.format(c).replace(/^−/, "-");
    let negative = false;
    const simple = (t) => !/[\s+\-−]/.test(t);
    if (text.startsWith("-") && simple(text.slice(1))) {
      negative = true;
      text = text.slice(1);
    }
    const isOne = text === "1";
    let coeff = text;
    if (!simple(text)) coeff = mono ? `(${text})` : text;
    else {
      const frac = /^(\d+)\/(\d+)$/.exec(text);
      if (frac && style === "tex") coeff = `\\frac{${frac[1]}}{${frac[2]}}`;
      else if (frac && style === "unicode" && mono) coeff = `(${text})`;
    }
    let body;
    if (!mono) body = coeff;
    else if (isOne) body = mono;
    else if (style === "text") body = `${coeff}*${mono}`;
    else body = `${coeff}${/^(\d+|\(.*\)|\\frac\{\d+\}\{\d+\})$/.test(coeff) ? "" : style === "tex" ? "\\," : "·"}${mono}`;
    if (idx === 0) out.push(negative ? `${minus}${body}` : body);
    else out.push(negative ? ` ${minus} ${body}` : ` + ${body}`);
  });
  return out.join("");
}
