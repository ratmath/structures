// Coefficient fields.
//
// A field is the plain object the rest of RatMath already uses
// (@ratmath/linalg's `field` module, @ratmath/algebraic's primeField,
// galoisField and NumberField):
//
//   { zero, one, add, sub, mul, div, neg, sign(v) -> 0 | 1 }
//
// Polynomial code only ever touches coefficients through such an object, so
// one implementation serves Q, F_p, GF(p^k) and Q(α). `normalize` fills in
// the few extras the polynomial code wants (from, equals, format, toPublic).

import { Rational } from "@ratmath/core";
import { toPair, addQ, subQ, mulQ, divQ, pairToString } from "@ratmath/algebraic/q";

const ZERO = [0n, 1n];
const ONE = [1n, 1n];

/**
 * The rationals. Elements are [numerator, denominator] BigInt pairs inside
 * the library (fast, no allocation of classes in inner loops) and Core
 * `Rational`s wherever a value is handed back to the caller.
 */
export const QQ = Object.freeze({
  name: "Q",
  characteristic: 0n,
  zero: ZERO,
  one: ONE,
  add: addQ,
  sub: subQ,
  mul: mulQ,
  div: (a, b) => {
    if (b[0] === 0n) throw new RangeError("division by zero");
    return divQ(a, b);
  },
  neg: ([n, d]) => [-n, d],
  inv: ([n, d]) => {
    if (n === 0n) throw new RangeError("division by zero");
    return n < 0n ? [-d, -n] : [d, n];
  },
  sign: ([n]) => (n === 0n ? 0 : 1),
  equals: (a, b) => a[0] === b[0] && a[1] === b[1],
  from: (x) => (x instanceof Rational ? [x.numerator, x.denominator] : toPair(x)),
  format: pairToString,
  toPublic: ([n, d]) => new Rational(n, d),
  /** −1, 0 or 1: Q is ordered, which printing uses to write a − b rather than a + −b. */
  order: ([n]) => (n < 0n ? -1 : n > 0n ? 1 : 0),
});

const formatAny = (v) => {
  if (Array.isArray(v)) return `(${v.map(String).join(", ")})`;
  return String(v);
};

const cache = new WeakMap();

function fieldName(field) {
  if (field.m && field.n && typeof field.name === "string") return `Q(${field.name})`; // a NumberField
  if (field.p !== undefined) return field.k > 1 ? `GF(${field.p}^${field.k})` : `F_${field.p}`;
  return typeof field.name === "string" ? field.name : "K";
}

/**
 * Fill in the optional members of a field object so the polynomial code can
 * rely on them. The object handed in is not changed.
 */
export function normalize(field) {
  if (field === QQ || field === undefined || field === null) return QQ;
  if (typeof field !== "object") throw new TypeError("a field must be an object with zero, one, add, sub, mul, div, neg and sign");
  const hit = cache.get(field);
  if (hit) return hit;
  for (const k of ["add", "sub", "mul", "div", "neg", "sign"]) {
    if (typeof field[k] !== "function") throw new TypeError(`this field has no ${k}()`);
  }
  if (!("zero" in field) || !("one" in field)) throw new TypeError("a field needs zero and one");
  const K = {
    name: fieldName(field),
    characteristic: field.characteristic ?? field.p ?? 0n,
    zero: field.zero,
    one: field.one,
    add: field.add,
    sub: field.sub,
    mul: field.mul,
    div: field.div,
    neg: field.neg,
    inv: field.inv ?? ((a) => field.div(field.one, a)),
    sign: (a) => {
      const s = field.sign(a);
      if (s === null || s === undefined) throw new RangeError("this field cannot decide whether a coefficient is zero; polynomials need exact coefficients");
      return s === 0 ? 0 : 1;
    },
    equals: field.equals ?? ((a, b) => field.sign(field.sub(a, b)) === 0),
    from: field.from ? (x) => field.from(liftInput(x)) : field.of ? (x) => field.of(liftInput(x)) : (x) => x,
    format: field.format ?? formatAny,
    toPublic: (a) => a,
    order: null,
    source: field,
  };
  cache.set(field, K);
  return K;
}

// Text such as "3/4" is read as an exact rational before a field's own
// `from` sees it, so every field accepts the same literals.
function liftInput(x) {
  if (typeof x === "string") {
    const [n, d] = toPair(x);
    return new Rational(n, d);
  }
  return x;
}

/** An integer n as an element of K (n·1), by doubling so that it works in any characteristic. */
export function fromInteger(K, n) {
  if (K === QQ) return [BigInt(n), 1n];
  n = BigInt(n);
  const neg = n < 0n;
  if (neg) n = -n;
  let acc = K.zero;
  let pow = K.one;
  while (n > 0n) {
    if (n & 1n) acc = K.add(acc, pow);
    pow = K.add(pow, pow);
    n >>= 1n;
  }
  return neg ? K.neg(acc) : acc;
}
