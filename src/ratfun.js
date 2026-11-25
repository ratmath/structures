// Rational functions: fractions of polynomials in several variables, always
// kept in lowest terms. Over Q both parts have integer coefficients with no
// common factor and the denominator leads with a positive coefficient; over
// other fields the denominator's leading coefficient is 1.
//
//   const R = ring("s, t");
//   const F = fractions(R);                 // the field Q(s, t)
//   F.from("(s^2 - t^2)/(s - t)")           // s + t
//
// `fractions(R)` is a field object in the shape the rest of RatMath uses
// ({ zero, one, add, sub, mul, div, neg, sign }), so @ratmath/linalg's
// field-generic elimination, determinant and inverse work on matrices whose
// entries contain parameters.

import { MPoly } from "./mpoly.js";
import { gcd } from "./gcd.js";
import { parseWith } from "./parse.js";
import { QQ } from "./field.js";

const gcdInt = (a, b) => {
  if (a < 0n) a = -a;
  if (b < 0n) b = -b;
  while (b) [a, b] = [b, a % b];
  return a;
};

export class RationalFunction {
  /** Use `fractions(ring).from(...)` or `.of(num, den)`; the constructor expects reduced parts. */
  constructor(field, num, den) {
    this.field = field;
    this.num = num;
    this.den = den;
    Object.freeze(this);
  }

  get ring() {
    return this.field.ring;
  }

  isZero() {
    return this.num.isZero();
  }

  /** Whether the denominator is a constant, so this is a polynomial. */
  isPolynomial() {
    return this.den.isConstant();
  }

  /** The polynomial this is; throws if the denominator is not a constant. */
  toPolynomial() {
    if (!this.isPolynomial()) throw new RangeError(`${this} is not a polynomial`);
    // over Q the parts are kept with integer coefficients, so (x + 1)/2 has denominator 2
    return this.den.equals(this.ring.one) ? this.num : this.num.div(this.den);
  }

  add(other) {
    const b = this.field.from(other);
    if (this.den.equals(b.den)) return this.field.of(this.num.add(b.num), this.den);
    // cross-multiply through the lcm of the denominators only
    const g = gcd(this.den, b.den);
    const d1 = this.den.div(g);
    const d2 = b.den.div(g);
    return this.field.of(this.num.mul(d2).add(b.num.mul(d1)), d1.mul(b.den));
  }

  neg() {
    return new RationalFunction(this.field, this.num.neg(), this.den);
  }

  sub(other) {
    return this.add(this.field.from(other).neg());
  }

  mul(other) {
    const b = this.field.from(other);
    if (this.isZero() || b.isZero()) return this.field.zero;
    // cancel across before multiplying: both fractions are already reduced
    const g1 = gcd(this.num, b.den);
    const g2 = gcd(b.num, this.den);
    return this.field.of(this.num.div(g1).mul(b.num.div(g2)), this.den.div(g2).mul(b.den.div(g1)), true);
  }

  inv() {
    if (this.isZero()) throw new RangeError("division by zero");
    return this.field.of(this.den, this.num, true);
  }

  div(other) {
    return this.mul(this.field.from(other).inv());
  }

  pow(n) {
    if (typeof n === "bigint") n = Number(n);
    if (!Number.isInteger(n)) throw new RangeError("a rational function can only be raised to a whole-number power");
    if (n < 0) return this.inv().pow(-n);
    return new RationalFunction(this.field, this.num.pow(n), this.den.pow(n));
  }

  equals(other) {
    const b = this.field.from(other);
    return this.num.equals(b.num) && this.den.equals(b.den);
  }

  /** The partial derivative with respect to a variable, by the quotient rule. */
  diff(v) {
    return this.field.of(this.num.diff(v).mul(this.den).sub(this.num.mul(this.den.diff(v))), this.den.mul(this.den));
  }

  /** Substitute polynomials or coefficients for some variables; throws if the denominator becomes zero. */
  subs(values) {
    const d = this.den.subs(values);
    if (d.isZero()) throw new RangeError(`the denominator of ${this} vanishes there`);
    return this.field.of(this.num.subs(values), d);
  }

  /** The value at a point (every variable that occurs must be given); throws at a pole. */
  at(point) {
    const K = this.ring.K;
    const d = K.from(this.den.at(point));
    if (K.sign(d) === 0) throw new RangeError(`${this} has a pole there`);
    return K.toPublic(K.div(K.from(this.num.at(point)), d));
  }

  toString() {
    return formatFraction(this, (p) => p.toString(), (a, b) => `${a}/${b}`);
  }

  toText() {
    return formatFraction(this, (p) => p.toText(), (a, b) => `${a}/${b}`);
  }

  toTeX() {
    if (this.isPolynomial()) return this.toPolynomial().toTeX();
    return `\\frac{${this.num.toTeX()}}{${this.den.toTeX()}}`;
  }
}

function formatFraction(f, show, join) {
  if (f.isPolynomial()) return show(f.toPolynomial());
  const wrap = (p) => (p.length > 1 ? `(${show(p)})` : show(p));
  const den = f.den.length > 1 || !f.den.isMonomial() || f.den.degree() > 1 || f.den.variables().length > 1 ? `(${show(f.den)})` : show(f.den);
  return join(wrap(f.num), den);
}

const fieldCache = new WeakMap();

/**
 * The field of fractions of a polynomial ring: rational functions in its
 * variables. The result is a field object usable wherever RatMath takes one.
 */
export function fractions(R) {
  const hit = fieldCache.get(R);
  if (hit) return hit;
  const F = {
    ring: R,
    name: `${R.K.name}(${R.vars.join(", ")})`,
    characteristic: R.K.characteristic,
    /** A fraction num/den, reduced to lowest terms. */
    of(num, den = R.one, reduced = false) {
      num = R.from(num);
      den = R.from(den);
      if (den.isZero()) throw new RangeError("division by zero");
      if (num.isZero()) return F.zero;
      if (!reduced) {
        const g = gcd(num, den);
        if (!g.isConstant()) {
          num = num.div(g);
          den = den.div(g);
        }
      }
      if (R.K === QQ) {
        // over Q: integer coefficients with no common factor, denominator leading positively
        let l = 1n;
        let g = 0n;
        for (const p of [num, den]) for (const { c } of p._m.values()) l = (l * c[1]) / gcdInt(l, c[1]);
        for (const p of [num, den]) for (const { c } of p._m.values()) g = gcdInt(g, (c[0] * l) / c[1]);
        const scale = den._lead().c[0] < 0n ? [-l, g] : [l, g];
        if (scale[0] !== scale[1]) {
          num = num._scale(scale);
          den = den._scale(scale);
        }
      } else {
        const lead = den._lead().c;
        if (!R.K.equals(lead, R.K.one)) {
          const inv = R.K.inv(lead);
          num = num._scale(inv);
          den = den._scale(inv);
        }
      }
      return new RationalFunction(F, num, den);
    },
    /** Anything reasonable as a rational function: one already, a polynomial, a coefficient, or text. */
    from(x) {
      if (x instanceof RationalFunction) {
        if (x.field === F) return x;
        return F.of(x.num.into(R), x.den.into(R));
      }
      if (x instanceof MPoly) return F.of(x);
      if (typeof x === "string") return F.parse(x);
      return F.of(R.from(x));
    },
    /** Read text such as "(x^2 - 1)/(x y + 1)"; negative exponents are allowed. */
    parse(text) {
      return parseWith(
        {
          vars: R.vars,
          number: (t) => F.of(R.constant(R.K.from(t))),
          variable: (v) => F.gens[R.indexOf(v)],
          integer: (f) => {
            if (!f.isPolynomial() || !f.num.isConstant()) return null;
            const c = f.isZero() ? "0" : R.K.format(f.num._lead().c);
            return /^-?\d+$/.test(c) ? Number(c) : null;
          },
          divide: (a, b) => a.div(b),
          inexact: "",
          negativePowers: true,
          name: F.name,
        },
        text,
      );
    },
    add: (a, b) => F.from(a).add(b),
    sub: (a, b) => F.from(a).sub(b),
    mul: (a, b) => F.from(a).mul(b),
    div: (a, b) => F.from(a).div(b),
    neg: (a) => F.from(a).neg(),
    inv: (a) => F.from(a).inv(),
    pow: (a, n) => F.from(a).pow(n),
    sign: (a) => (F.from(a).isZero() ? 0 : 1),
    equals: (a, b) => F.from(a).equals(b),
    format: (a) => F.from(a).toString(),
    toString: () => F.name,
  };
  F.zero = new RationalFunction(F, R.zero, R.one);
  F.one = new RationalFunction(F, R.one, R.one);
  F.gens = Object.freeze(R.gens.map((g) => new RationalFunction(F, g, R.one)));
  fieldCache.set(R, Object.freeze(F));
  return F;
}
