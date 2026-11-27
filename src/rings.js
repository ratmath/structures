// The rings that are met first: the integers modulo n, and the rings of
// integers of quadratic fields (the Gaussian and Eisenstein integers among
// them), with their units, Euclidean division where the norm allows it,
// factorisation of elements, and the ideals that put unique factorisation
// back where elements lose it.
//
//   const Z12 = Zn(12);
//   Z12.units();  Z12.zeroDivisors();  Z12.idempotents();  Z12.decompose();
//   const R = quadraticRing(-5);                 // Z[√−5]
//   R.factorizations(6);                         // 2·3 and (1 + √−5)(1 − √−5)
//   R.ideal(2, R.of(1, 1)).isPrincipal();        // false
//   R.classNumber();                             // 2

import { factorInteger, totient, modInverse, modPow, isProbablePrime, sqrtMod, primeField } from "@ratmath/algebraic";

const B = (x) => {
  if (typeof x === "bigint") return x;
  if (typeof x === "number") {
    if (!Number.isInteger(x)) throw new RangeError(`${x} is not a whole number`);
    return BigInt(x);
  }
  if (typeof x === "string" && /^-?\d+$/.test(x.trim())) return BigInt(x.trim());
  if (x && typeof x === "object" && "numerator" in x && x.denominator === 1n) return x.numerator;
  throw new RangeError(`${String(x)} is not a whole number`);
};
const abs = (a) => (a < 0n ? -a : a);
const gcd = (a, b) => {
  a = abs(a);
  b = abs(b);
  while (b) [a, b] = [b, a % b];
  return a;
};
const mod = (a, n) => ((a % n) + n) % n;
function egcd(a, b) {
  let [r0, r1, s0, s1, t0, t1] = [a, b, 1n, 0n, 0n, 1n];
  while (r1) {
    const q = r0 / r1;
    [r0, r1] = [r1, r0 - q * r1];
    [s0, s1] = [s1, s0 - q * s1];
    [t0, t1] = [t1, t0 - q * t1];
  }
  return r0 < 0n ? [-r0, -s0, -t0] : [r0, s0, t0];
}
function isqrt(n) {
  if (n < 0n) return null;
  if (n < 2n) return n;
  let x = BigInt(Math.floor(Math.sqrt(Number(n))));
  while (x * x > n) x -= 1n;
  while ((x + 1n) * (x + 1n) <= n) x += 1n;
  return x;
}
const exactSqrt = (n) => {
  const r = isqrt(n);
  return r !== null && r * r === n ? r : null;
};
const divisors = (n) => {
  let out = [1n];
  for (const { prime, exponent } of factorInteger(n)) {
    const next = [];
    let pk = 1n;
    for (let k = 0; k <= exponent; k += 1) {
      for (const d of out) next.push(d * pk);
      pk *= prime;
    }
    out = next;
  }
  return out.sort((a, b) => (a < b ? -1 : a > b ? 1 : 0));
};

/**
 * The Chinese remainder theorem for integers: the x with x ≡ rᵢ (mod mᵢ)
 * for every i. The moduli need not be coprime; when the congruences
 * contradict each other the result is null. Returns { value, modulus }.
 */
export function crt(residues, moduli) {
  if (residues.length !== moduli.length || !residues.length) throw new RangeError("give as many residues as moduli, at least one");
  let x = 0n;
  let m = 1n;
  for (let i = 0; i < residues.length; i += 1) {
    const mi = abs(B(moduli[i]));
    if (mi === 0n) throw new RangeError("a modulus is zero");
    const ri = mod(B(residues[i]), mi);
    const [g, s] = egcd(m, mi);
    if ((ri - x) % g !== 0n) return null;
    const l = (m / g) * mi;
    x = mod(x + m * mod(((ri - x) / g) * s, mi / g), l);
    m = l;
  }
  return { value: x, modulus: m };
}

// ---- Z/n -----------------------------------------------------------------------------------

/** The ring of integers modulo n. Elements are the BigInts 0, …, n − 1. */
export class ZnRing {
  constructor(n) {
    this.n = B(n);
    if (this.n < 1n) throw new RangeError("the modulus must be a positive whole number");
    this.factorization = factorInteger(this.n);
    this.zero = 0n;
    this.one = mod(1n, this.n);
  }
  get characteristic() {
    return this.n;
  }
  get size() {
    return this.n;
  }
  of(a) {
    return mod(B(a), this.n);
  }
  add(a, b) {
    return mod(B(a) + B(b), this.n);
  }
  sub(a, b) {
    return mod(B(a) - B(b), this.n);
  }
  mul(a, b) {
    return mod(B(a) * B(b), this.n);
  }
  neg(a) {
    return mod(-B(a), this.n);
  }
  pow(a, e) {
    e = B(e);
    return e < 0n ? modPow(this.inverse(a), -e, this.n) : modPow(this.of(a), e, this.n);
  }
  isUnit(a) {
    return gcd(this.of(a), this.n) === 1n;
  }
  inverse(a) {
    if (!this.isUnit(a)) throw new RangeError(`${this.of(a)} has no inverse modulo ${this.n}: it shares the factor ${gcd(this.of(a), this.n)} with the modulus`);
    return this.n === 1n ? 0n : modInverse(this.of(a), this.n);
  }
  div(a, b) {
    return this.mul(a, this.inverse(b));
  }
  isZeroDivisor(a) {
    const x = this.of(a);
    return this.n > 1n && gcd(x, this.n) !== 1n;
  }
  /** a field exactly when n is prime */
  get isField() {
    return isProbablePrime(this.n);
  }
  get isIntegralDomain() {
    return this.isField;
  }
  /** no nilpotent elements: n is square-free */
  get isReduced() {
    return this.factorization.every((f) => f.exponent === 1);
  }
  _small(what) {
    if (this.n > 100000n) throw new RangeError(`listing ${what} needs a modulus of at most 100000`);
  }
  elements() {
    this._small("the elements");
    return Array.from({ length: Number(this.n) }, (_, i) => BigInt(i));
  }
  units() {
    return this.elements().filter((a) => this.isUnit(a));
  }
  /** Euler's φ(n): the number of units */
  get unitCount() {
    return this.n === 1n ? 1n : totient(this.n);
  }
  zeroDivisors() {
    return this.elements().filter((a) => a !== 0n && this.isZeroDivisor(a));
  }
  /** the elements with some power zero: the multiples of the product of the primes dividing n */
  nilpotents() {
    const rad = this.factorization.reduce((p, f) => p * f.prime, 1n);
    this._small("the nilpotent elements");
    return this.elements().filter((a) => a % rad === 0n);
  }
  /** the elements with e² = e: one for each way of splitting the prime powers of n into two groups */
  idempotents() {
    const parts = this.factorization.map((f) => f.prime ** BigInt(f.exponent));
    const out = [];
    for (let mask = 0; mask < 1 << parts.length; mask += 1) {
      out.push(crt(parts.map((_, i) => (mask & (1 << i) ? 1n : 0n)), parts).value % this.n);
    }
    return out.sort((a, b) => (a < b ? -1 : 1));
  }
  /** the order of a in the additive group */
  additiveOrder(a) {
    return this.n / gcd(this.of(a), this.n);
  }
  /** the order of a unit in the multiplicative group */
  orderOf(a) {
    const x = this.of(a);
    if (!this.isUnit(x)) throw new RangeError(`${x} is not a unit modulo ${this.n}`);
    if (this.n === 1n) return 1n;
    let ord = this.unitCount;
    for (const { prime } of factorInteger(ord)) while (ord % prime === 0n && modPow(x, ord / prime, this.n) === 1n) ord /= prime;
    return ord;
  }
  /**
   * The group of units as a product of cyclic groups: the orders of the
   * factors, each dividing the next. One factor means the group is cyclic
   * (there is a primitive root): exactly for n = 1, 2, 4, pᵏ and 2pᵏ.
   */
  unitGroup() {
    const cyclic = [];
    for (const { prime: p, exponent: k } of this.factorization) {
      if (p === 2n) {
        if (k === 2) cyclic.push(2n);
        else if (k > 2) cyclic.push(2n, 2n ** BigInt(k - 2));
      } else cyclic.push((p - 1n) * p ** BigInt(k - 1));
    }
    // invariant factors from the prime-power parts
    const byPrime = new Map();
    for (const c of cyclic) for (const { prime, exponent } of factorInteger(c)) {
      if (!byPrime.has(prime)) byPrime.set(prime, []);
      byPrime.get(prime).push(prime ** BigInt(exponent));
    }
    const rank = Math.max(0, ...[...byPrime.values()].map((l) => l.length));
    const invariants = Array.from({ length: rank }, () => 1n);
    for (const list of byPrime.values()) {
      list.sort((a, b) => (a < b ? 1 : -1));
      list.forEach((q, i) => {
        invariants[rank - 1 - i] *= q;
      });
    }
    return { order: this.unitCount, invariants, cyclic: rank <= 1, text: invariants.length ? invariants.map((q) => `C${q}`).join(" × ") : "trivial" };
  }
  /** a generator of the group of units, or null when the group is not cyclic */
  primitiveRoot() {
    if (!this.unitGroup().cyclic) return null;
    const phi = this.unitCount;
    for (let a = 1n; a < this.n || this.n === 1n; a += 1n) {
      if (this.n === 1n) return 0n;
      if (this.isUnit(a) && this.orderOf(a) === phi) return a;
    }
    return null;
  }
  /** The ideals: one for each divisor d of n, the multiples of d. */
  ideals() {
    return divisors(this.n).map((d) => ({ generator: d % this.n, size: this.n / d, quotient: `Z/${d}`, isMaximal: isProbablePrime(d), isPrime: isProbablePrime(d) }));
  }
  /**
   * The Chinese remainder theorem: Z/n as the product of the rings Z/pᵏ for
   * the prime powers of n. `to` and `from` are the isomorphism and its
   * inverse, and `idempotents[i]` is the element that is 1 in factor i and
   * 0 in the others.
   */
  decompose() {
    const moduli = this.factorization.map((f) => f.prime ** BigInt(f.exponent));
    return {
      moduli,
      to: (a) => moduli.map((m) => mod(B(a), m)),
      from: (rs) => (moduli.length ? crt(rs, moduli).value : 0n),
      idempotents: moduli.map((_, i) => crt(moduli.map((__, j) => (i === j ? 1n : 0n)), moduli).value),
    };
  }
  /** the addition or multiplication table, for a small modulus */
  table(op = "mul") {
    if (this.n > 64n) throw new RangeError("a table needs a modulus of at most 64");
    const els = this.elements();
    return els.map((a) => els.map((b) => (op === "add" ? this.add(a, b) : this.mul(a, b))));
  }
  /** this ring as a field object (for ring(…, { field }) and linalg), when n is prime */
  asField() {
    if (!this.isField) throw new RangeError(`Z/${this.n} is not a field: ${this.n} is not prime`);
    return primeField(this.n);
  }
  toString() {
    return `Z/${this.n}`;
  }
}
export const Zn = (n) => new ZnRing(n);

// ---- quadratic integer rings ---------------------------------------------------------------

const NORM_EUCLIDEAN = [-11, -7, -3, -2, -1, 2, 3, 5, 6, 7, 11, 13, 17, 19, 21, 29, 33, 37, 41, 57, 73];

/** An element a + bω of a quadratic ring. */
export class QuadraticInteger {
  constructor(R, a, b) {
    this.ring = R;
    this.a = a;
    this.b = b;
  }
  /** twice the coordinates on 1 and √d: the element is (x + y√d)/2 */
  get xy() {
    return this.ring.half ? [2n * this.a + this.b, this.b] : [2n * this.a, 2n * this.b];
  }
  add(o) {
    o = this.ring.of(o);
    return new QuadraticInteger(this.ring, this.a + o.a, this.b + o.b);
  }
  sub(o) {
    o = this.ring.of(o);
    return new QuadraticInteger(this.ring, this.a - o.a, this.b - o.b);
  }
  neg() {
    return new QuadraticInteger(this.ring, -this.a, -this.b);
  }
  mul(o) {
    o = this.ring.of(o);
    // ω² = t·ω − n
    const { t, n } = this.ring;
    const bb = this.b * o.b;
    return new QuadraticInteger(this.ring, this.a * o.a - n * bb, this.a * o.b + this.b * o.a + t * bb);
  }
  pow(k) {
    let out = this.ring.one;
    let base = this;
    if (k < 0) {
      if (!this.isUnit()) throw new RangeError(`${this} is not a unit`);
      base = this.inverse();
      k = -k;
    }
    for (; k > 0; k >>= 1) {
      if (k & 1) out = out.mul(base);
      base = base.mul(base);
    }
    return out;
  }
  /** the conjugate: √d ↦ −√d */
  conj() {
    return new QuadraticInteger(this.ring, this.a + this.ring.t * this.b, -this.b);
  }
  /** the norm a² + t·ab + n·b²: the product with the conjugate, an integer */
  norm() {
    const { t, n } = this.ring;
    return this.a * this.a + t * this.a * this.b + n * this.b * this.b;
  }
  trace() {
    return 2n * this.a + this.ring.t * this.b;
  }
  isZero() {
    return this.a === 0n && this.b === 0n;
  }
  equals(o) {
    o = this.ring.of(o);
    return this.a === o.a && this.b === o.b;
  }
  isUnit() {
    return abs(this.norm()) === 1n;
  }
  inverse() {
    const N = this.norm();
    if (abs(N) !== 1n) throw new RangeError(`${this} is not a unit: its norm is ${N}`);
    const c = this.conj();
    return N === 1n ? c : c.neg();
  }
  /** the exact quotient, or null when o does not divide this */
  divExact(o) {
    o = this.ring.of(o);
    if (o.isZero()) throw new RangeError("division by zero");
    const N = o.norm();
    const p = this.mul(o.conj());
    return p.a % N === 0n && p.b % N === 0n ? new QuadraticInteger(this.ring, p.a / N, p.b / N) : null;
  }
  isDivisibleBy(o) {
    return this.divExact(o) !== null;
  }
  isAssociateOf(o) {
    o = this.ring.of(o);
    return !o.isZero() && this.isDivisibleBy(o) && o.isDivisibleBy(this);
  }
  /** −1, 0 or 1: the sign of the real number, for a real quadratic ring */
  sign() {
    if (this.ring.d < 0n) throw new RangeError("only elements of a real quadratic ring have a sign");
    const [x, y] = this.xy;
    const sx = x > 0n ? 1 : x < 0n ? -1 : 0;
    const sy = y > 0n ? 1 : y < 0n ? -1 : 0;
    if (sx === sy || sy === 0) return sx;
    if (sx === 0) return sy;
    const c = x * x - this.ring.d * y * y; // |x| against |y|√d
    return c > 0n ? sx : sy;
  }
  toNumber() {
    const [x, y] = this.xy;
    if (this.ring.d < 0n) throw new RangeError("not a real number");
    return (Number(x) + Number(y) * Math.sqrt(Number(this.ring.d))) / 2;
  }
  toString() {
    const R = this.ring;
    let [x, y] = this.xy;
    const halves = x % 2n !== 0n || y % 2n !== 0n;
    if (!halves) {
      x /= 2n;
      y /= 2n;
    }
    const root = R.d === -1n ? "i" : `√${String(R.d).replace("-", "−")}`;
    const num = (v) => String(v).replace("-", "−");
    let body;
    if (y === 0n) body = num(x);
    else {
      const ay = abs(y);
      const term = `${ay === 1n ? "" : ay}${root}`;
      body = x === 0n ? `${y < 0n ? "−" : ""}${term}` : `${num(x)} ${y < 0n ? "−" : "+"} ${term}`;
    }
    return halves ? `(${body})/2` : body;
  }
}

/** A non-zero ideal of a quadratic ring, as the lattice Z·a + Z·(b + cω). */
export class QuadraticIdeal {
  constructor(R, a, b, c) {
    this.ring = R;
    this.a = a;
    this.b = b;
    this.c = c;
  }
  /** the number of residue classes: the size of the quotient ring */
  get norm() {
    return this.a * this.c;
  }
  /** two elements that generate the ideal */
  get generators() {
    return [this.ring.of(this.a), new QuadraticInteger(this.ring, this.b, this.c)];
  }
  contains(x) {
    x = this.ring.of(x);
    if (x.b % this.c !== 0n) return false;
    return (x.a - (x.b / this.c) * this.b) % this.a === 0n;
  }
  equals(o) {
    return this.a === o.a && this.b === o.b && this.c === o.c;
  }
  isSubsetOf(o) {
    return this.generators.every((g) => o.contains(g));
  }
  add(o) {
    return this.ring.ideal(...this.generators, ...o.generators);
  }
  mul(o) {
    o = o instanceof QuadraticIdeal ? o : this.ring.ideal(o);
    const out = [];
    for (const g of this.generators) for (const h of o.generators) out.push(g.mul(h));
    return this.ring.ideal(...out);
  }
  pow(k) {
    let out = this.ring.ideal(1);
    for (let i = 0; i < k; i += 1) out = out.mul(this);
    return out;
  }
  conj() {
    return this.ring.ideal(...this.generators.map((g) => g.conj()));
  }
  get isWhole() {
    return this.norm === 1n;
  }
  /** a generator when the ideal is principal, otherwise null */
  generator() {
    if (this._gen !== undefined) return this._gen;
    for (const x of this.ring.elementsOfNorm(this.norm, { bothSigns: true })) {
      // every unit multiple of x has the same membership
      if (this.contains(x)) return (this._gen = x);
    }
    return (this._gen = null);
  }
  isPrincipal() {
    return this.generator() !== null;
  }
  /** prime exactly when the quotient is a field: the norm is a prime, or the ideal is (p) for an inert prime */
  isPrime() {
    const N = this.norm;
    if (isProbablePrime(N)) return true;
    const p = exactSqrt(N);
    return p !== null && isProbablePrime(p) && this.equals(this.ring.ideal(p)) && this.ring.splitting(p).kind === "inert";
  }
  /** is this ideal in the same class as another (I = α·J for some α of the field)? */
  isEquivalentTo(o) {
    return this.mul(o.conj()).isPrincipal();
  }
  /** The factorisation into prime ideals, which always exists and is unique: [{ prime, exponent }]. */
  factor() {
    const out = [];
    for (const { prime: p } of factorInteger(this.norm)) {
      for (const P of this.ring.primesAbove(p)) {
        let k = 0;
        let Q = P;
        while (this.isSubsetOf(Q)) {
          k += 1;
          Q = Q.mul(P);
        }
        if (k) out.push({ prime: P, exponent: k });
      }
    }
    return out;
  }
  toString() {
    const g = this.generator?.call(this);
    if (g) return `(${g})`;
    const [x, y] = this.generators;
    return `(${x}, ${y})`;
  }
}

/** The ring of integers of Q(√d), for a square-free integer d other than 0 and 1. */
export class QuadraticRing {
  constructor(d) {
    this.d = B(d);
    if (this.d === 0n || this.d === 1n) throw new RangeError("d must not be 0 or 1");
    for (const { exponent } of factorInteger(this.d)) if (exponent > 1) throw new RangeError(`${this.d} is not square-free; use its square-free part`);
    /** true when d ≡ 1 (mod 4): then ω = (1 + √d)/2, otherwise ω = √d */
    this.half = mod(this.d, 4n) === 1n;
    // ω² = t·ω − n
    this.t = this.half ? 1n : 0n;
    this.n = this.half ? (1n - this.d) / 4n : -this.d;
    this.discriminant = this.half ? this.d : 4n * this.d;
    this.zero = new QuadraticInteger(this, 0n, 0n);
    this.one = new QuadraticInteger(this, 1n, 0n);
    this.omega = new QuadraticInteger(this, 0n, 1n);
    /** √d as an element */
    this.sqrtD = this.half ? new QuadraticInteger(this, -1n, 2n) : this.omega;
    this.name = this.d === -1n ? "Z[i]" : this.half ? `Z[(1 + √${String(this.d).replace("-", "−")})/2]` : `Z[√${String(this.d).replace("-", "−")}]`;
  }
  /** An element: of(a, b) is a + b·ω; of(n) an integer; of("3 + 2s") is not read, use sqrt(x, y). */
  of(a, b) {
    if (a instanceof QuadraticInteger) {
      if (a.ring !== this) throw new TypeError("that element belongs to another ring");
      return a;
    }
    return new QuadraticInteger(this, B(a), b === undefined ? 0n : B(b));
  }
  /** x + y√d, for whole numbers x and y (or both halves of odd numbers when d ≡ 1 mod 4: sqrt("1/2", "1/2")) */
  sqrt(x, y) {
    const twice = (v) => {
      if (typeof v === "string" && /^-?\d+\/2$/.test(v.trim())) return BigInt(v.trim().split("/")[0]);
      if (v && typeof v === "object" && "numerator" in v && (v.denominator === 1n || v.denominator === 2n)) return v.denominator === 1n ? 2n * v.numerator : v.numerator;
      return 2n * B(v);
    };
    const [X, Y] = [twice(x), twice(y)];
    // (X + Y√d)/2 = a + bω
    if (this.half) {
      if (mod(X - Y, 2n) !== 0n) throw new RangeError("not an integer of this ring: x and y must both be whole, or both halves of odd numbers");
      return new QuadraticInteger(this, (X - Y) / 2n, Y);
    }
    if (X % 2n !== 0n || Y % 2n !== 0n) throw new RangeError("not an integer of this ring: x and y must be whole numbers");
    return new QuadraticInteger(this, X / 2n, Y / 2n);
  }
  get isReal() {
    return this.d > 0n;
  }
  /**
   * The units. Imaginary rings have finitely many (4 in Z[i], 6 in the
   * Eisenstein integers, otherwise ±1). Real rings have ±εᵏ for a
   * fundamental unit ε > 1, found from the continued fraction of ω.
   */
  get units() {
    if (this.d === -1n) return { finite: true, list: [this.one, this.omega, this.one.neg(), this.omega.neg()], count: 4 };
    if (this.d === -3n) {
      const w = this.omega;
      const list = [0, 1, 2, 3, 4, 5].map((k) => w.pow(k));
      return { finite: true, list, count: 6 };
    }
    if (this.d < 0n) return { finite: true, list: [this.one, this.one.neg()], count: 2 };
    const e = this.fundamentalUnit;
    return { finite: false, fundamental: e, norm: e.norm(), text: `±(${e})^k` };
  }
  /** the smallest unit greater than 1 of a real quadratic ring */
  get fundamentalUnit() {
    if (this.d < 0n) throw new RangeError("an imaginary quadratic ring has finitely many units");
    if (this._unit) return this._unit;
    // continued fraction of ω = (P + √d)/Q
    const D = this.d;
    const s = isqrt(D);
    let [P, Qd] = this.half ? [1n, 2n] : [0n, 1n];
    let [p0, p1, q0, q1] = [0n, 1n, 1n, 0n];
    for (let k = 0; k < 2000000; k += 1) {
      const a = (P + s) / Qd;
      [p0, p1] = [p1, a * p1 + p0];
      [q0, q1] = [q1, a * q1 + q0];
      // the convergent p/q approximates ω, so p − q·ω′ is large: a candidate unit
      const u = this.half ? new QuadraticInteger(this, p1 - q1, q1) : new QuadraticInteger(this, p1, q1);
      if (abs(u.norm()) === 1n) return (this._unit = u);
      P = a * Qd - P;
      Qd = (D - P * P) / Qd;
    }
    throw new RangeError("the fundamental unit was not found within the step limit");
  }
  /** Is the ring Euclidean for the absolute value of the norm? (The complete list of such d is known.) */
  get isNormEuclidean() {
    return NORM_EUCLIDEAN.includes(Number(this.d));
  }
  /**
   * Division with remainder: a = q·b + r with q the lattice point nearest
   * a/b. `euclidean` says whether |N(r)| < |N(b)|, which always holds in a
   * norm-Euclidean ring and can fail elsewhere.
   */
  divmod(a, b) {
    a = this.of(a);
    b = this.of(b);
    if (b.isZero()) throw new RangeError("division by zero");
    const N = b.norm();
    const p = a.mul(b.conj());
    // p/N in coordinates on 1 and ω; round each, then look at the neighbours
    const round = (x) => {
      const q = (2n * x + N) / (2n * N);
      return (2n * x + N) % (2n * N) < 0n ? q - 1n : q;
    };
    const [a0, b0] = [round(p.a), round(p.b)];
    let best = null;
    for (const da of [-1n, 0n, 1n]) for (const db of [-1n, 0n, 1n]) {
      const q = new QuadraticInteger(this, a0 + da, b0 + db);
      const r = a.sub(q.mul(b));
      const size = abs(r.norm());
      if (!best || size < best.size) best = { quotient: q, remainder: r, size };
    }
    return { quotient: best.quotient, remainder: best.remainder, euclidean: best.size < abs(N) };
  }
  /**
   * The greatest common divisor by Euclid's algorithm, with the steps and
   * Bézout coefficients: g = s·a + t·b. Refused when a division fails to
   * make the norm smaller, which happens only outside the norm-Euclidean rings.
   */
  gcd(a, b) {
    a = this.of(a);
    b = this.of(b);
    let [r0, r1, s0, s1, t0, t1] = [a, b, this.one, this.zero, this.zero, this.one];
    const steps = [];
    while (!r1.isZero()) {
      const { quotient, remainder, euclidean } = this.divmod(r0, r1);
      if (!euclidean) throw new RangeError(`Euclid's algorithm stalls in ${this.name}: dividing ${r0} by ${r1} leaves a remainder whose norm is not smaller. The ring is not Euclidean for the norm; use ideals instead`);
      steps.push({ dividend: r0, divisor: r1, quotient, remainder });
      [r0, r1] = [r1, remainder];
      [s0, s1] = [s1, s0.sub(quotient.mul(s1))];
      [t0, t1] = [t1, t0.sub(quotient.mul(t1))];
    }
    if (r0.isZero()) return { gcd: r0, s: s0, t: t0, steps };
    // hand back the preferred associate, and move the unit into s and t
    const g = this.canonical(r0);
    const u = g.divExact(r0);
    return { gcd: g, s: s0.mul(u), t: t0.mul(u), steps };
  }
  /** the comparison of two real elements as real numbers */
  _cmp(x, y) {
    return x.sub(y).sign();
  }
  /** the associate of x chosen as the representative of all its unit multiples */
  canonical(x) {
    x = this.of(x);
    if (x.isZero()) return x;
    if (this.d < 0n) {
      let best = null;
      for (const u of this.units.list) {
        const y = x.mul(u);
        if (!best || y.a > best.a || (y.a === best.a && y.b > best.b)) best = y;
      }
      return best;
    }
    // real: make it positive, then move by powers of ε until |m| ≤ y² < |m|·ε²
    const e = this.fundamentalUnit;
    const einv = e.inverse();
    let y = x.sign() < 0 ? x.neg() : x;
    const m = this.of(abs(x.norm()));
    const e2m = e.mul(e).mul(m);
    const guess = Math.floor((Math.log(Math.abs(Number(m.a)) || 1) / 2 - Math.log(Math.abs(y.toNumber()) || 1e-300)) / Math.log(e.toNumber()));
    if (Number.isFinite(guess) && guess !== 0) {
      y = y.mul(guess > 0 ? e.pow(guess) : einv.pow(-guess));
      if (y.sign() < 0) y = y.neg();
    }
    while (this._cmp(y.mul(y), m) < 0) {
      y = y.mul(e);
      if (y.sign() < 0) y = y.neg();
    }
    while (this._cmp(y.mul(y), e2m) >= 0) {
      y = y.mul(einv);
      if (y.sign() < 0) y = y.neg();
    }
    return y;
  }
  /**
   * The elements of norm m, one from each class of associates. With
   * { bothSigns: true }, norm −m as well. The search is finite: an
   * imaginary ring has finitely many elements of a given norm, and in a
   * real ring every class has a member of bounded size.
   */
  elementsOfNorm(m, { bothSigns = false, limit = 20000000 } = {}) {
    m = B(m);
    const d = this.d;
    const targets = bothSigns && d > 0n && m !== 0n ? [m, -m] : [m];
    const found = new Map();
    const push = (x) => {
      const c = this.canonical(x);
      // the preferred associate can have norm −m when there is a unit of norm −1: move back by that unit
      found.set(`${c.a},${c.b}`, !bothSigns && c.norm() !== m ? c.mul(this.fundamentalUnit) : c);
    };
    const am = abs(m);
    // (X + Y√d)/2 with X² − dY² = 4m
    let bound;
    if (d < 0n) bound = isqrt((4n * am) / -d);
    else {
      const e = this.fundamentalUnit.toNumber();
      const est = (2 * Math.sqrt(Number(am)) * (e + 1)) / Math.sqrt(Number(d)) + 2;
      if (!(est < limit)) throw new RangeError(`the search for elements of norm ${m} in ${this.name} is too long (the fundamental unit is about ${e.toExponential(2)})`);
      bound = BigInt(Math.ceil(est));
    }
    for (let Y = 0n; Y <= bound; Y += 1n) {
      for (const t of targets) {
        const X = exactSqrt(4n * t + d * Y * Y);
        if (X === null) continue;
        for (const sx of X === 0n ? [0n] : [X, -X]) {
          if (this.half ? mod(sx - Y, 2n) !== 0n : sx % 2n !== 0n || Y % 2n !== 0n) continue;
          const x = this.half ? new QuadraticInteger(this, (sx - Y) / 2n, Y) : new QuadraticInteger(this, sx / 2n, Y / 2n);
          if (!x.isZero() || m === 0n) push(x);
        }
      }
    }
    return [...found.values()];
  }
  /** Is x irreducible: not a unit, and not a product of two non-units? */
  isIrreducible(x) {
    x = this.of(x);
    const N = abs(x.norm());
    if (N <= 1n) return false;
    for (const m of divisors(N)) {
      if (m === 1n || m === N || m * m > N) continue;
      for (const y of this.elementsOfNorm(m, { bothSigns: true })) if (x.isDivisibleBy(y)) return false;
    }
    return true;
  }
  /** Is x prime: does it divide a product only when it divides a factor? Exactly when the ideal (x) is prime. */
  isPrimeElement(x) {
    x = this.of(x);
    if (x.isZero() || x.isUnit()) return false;
    return this.ideal(x).isPrime();
  }
  /** the irreducible divisors of x, one from each class of associates */
  irreducibleDivisors(x) {
    x = this.of(x);
    const N = abs(x.norm());
    const out = [];
    for (const m of divisors(N)) {
      if (m === 1n) continue;
      for (const y of this.elementsOfNorm(m, { bothSigns: true })) if (x.isDivisibleBy(y) && this.isIrreducible(y)) out.push(y);
    }
    return out;
  }
  /**
   * Every factorisation of x into irreducibles, up to order and units:
   * [{ unit, factors }] with x = unit·Π factors. One entry means the
   * factorisation of this element is unique; in a ring of class number
   * above 1 some elements have several.
   */
  factorizations(x, { limit = 200 } = {}) {
    x = this.of(x);
    if (x.isZero()) throw new RangeError("zero has no factorisation");
    const key = (y) => `${y.a},${y.b}`;
    const seen = new Set();
    const out = [];
    const irr = this.irreducibleDivisors(x);
    const walk = (rest, from, acc) => {
      if (out.length >= limit) return;
      if (rest.isUnit()) {
        const k = acc.map(key).sort().join(";");
        if (!seen.has(k)) {
          seen.add(k);
          out.push({ unit: rest, factors: [...acc] });
        }
        return;
      }
      for (let i = from; i < irr.length; i += 1) {
        const q = rest.divExact(irr[i]);
        if (q) walk(q, i, [...acc, irr[i]]);
      }
    };
    walk(x, 0, []);
    return out;
  }
  /** The ideal generated by some elements (or integers). */
  ideal(...gens) {
    const els = gens.flat().map((g) => this.of(g)).filter((g) => !g.isZero());
    if (!els.length) throw new RangeError("the zero ideal is not represented; give a non-zero generator");
    // the lattice spanned by g and g·ω for each generator, in Hermite normal form
    let a = 0n; // generator of the part with no ω
    let v = null; // (b, c) with c > 0 least
    const feed = (x, y) => {
      if (y === 0n) {
        a = gcd(a, x);
        return;
      }
      if (!v) {
        v = y < 0n ? [-x, -y] : [x, y];
        return;
      }
      const [g, s, t] = egcd(v[1], y);
      const kernel = (v[1] / g) * x - (y / g) * v[0];
      v = [s * v[0] + t * x, g];
      a = gcd(a, kernel);
    };
    for (const g of els) {
      feed(g.a, g.b);
      const w = g.mul(this.omega);
      feed(w.a, w.b);
    }
    return new QuadraticIdeal(this, a, mod(v[0], a), v[1]);
  }
  /**
   * How a rational prime p behaves: "inert" (it stays prime), "split" (it
   * is the product of two different prime ideals) or "ramified" (the square
   * of one), read from the roots of ω's polynomial modulo p.
   */
  splitting(p) {
    p = B(p);
    if (!isProbablePrime(p)) throw new RangeError(`${p} is not a prime`);
    const roots = this._roots(p);
    const kind = roots.length === 0 ? "inert" : roots.length === 1 ? "ramified" : "split";
    const primes = kind === "inert" ? [this.ideal(p)] : roots.map((r) => this.ideal(p, this.omega.sub(r)));
    return { prime: p, kind, primes, principal: primes.every((P) => P.isPrincipal()) };
  }
  /** the roots of x² − t·x + n modulo p */
  _roots(p) {
    const f = (x) => mod(x * x - this.t * x + this.n, p);
    if (p < 200n) {
      const out = [];
      for (let x = 0n; x < p; x += 1n) if (f(x) === 0n) out.push(x);
      return out;
    }
    const D = mod(this.discriminant, p);
    if (D === 0n) return [mod(this.t * modInverse(2n, p), p)];
    const s = sqrtMod(D, p);
    if (s === null || s === undefined) return [];
    const half = modInverse(2n, p);
    const rs = [mod((this.t + s) * half, p), mod((this.t - s) * half, p)];
    return rs.filter((r) => f(r) === 0n).sort((x, y) => (x < y ? -1 : 1));
  }
  /** the prime ideals containing p */
  primesAbove(p) {
    return this.splitting(p).primes;
  }
  /** Minkowski's bound: every ideal class contains an ideal of norm at most this */
  get minkowskiBound() {
    const D = Number(abs(this.discriminant));
    return this.d > 0n ? Math.sqrt(D) / 2 : (2 / Math.PI) * Math.sqrt(D);
  }
  /**
   * The class group, by Minkowski's bound: representatives of the ideal
   * classes among the ideals of small norm, and the class number h. h = 1
   * exactly when every ideal is principal, which here is the same as
   * unique factorisation of elements.
   */
  classGroup() {
    if (this._classes) return this._classes;
    const M = Math.floor(this.minkowskiBound);
    const primes = [];
    for (let p = 2; p <= M; p += 1) if (isProbablePrime(BigInt(p))) for (const P of this.primesAbove(p)) if (P.norm <= BigInt(M)) primes.push(P);
    // all ideals of norm ≤ M, as products of those primes
    const ideals = [this.ideal(1)];
    const seen = new Set(["1,0,1"]);
    const grow = (I, from) => {
      for (let i = from; i < primes.length; i += 1) {
        const J = I.mul(primes[i]);
        if (J.norm > BigInt(M)) continue;
        const k = `${J.a},${J.b},${J.c}`;
        if (seen.has(k)) continue;
        seen.add(k);
        ideals.push(J);
        grow(J, i);
      }
    };
    grow(ideals[0], 0);
    const classes = [];
    for (const I of ideals.sort((x, y) => (x.norm < y.norm ? -1 : x.norm > y.norm ? 1 : 0))) {
      if (!classes.some((c) => c.isEquivalentTo(I))) classes.push(I);
    }
    // the shape of the group, read off from the orders of the classes
    const h = classes.length;
    const orders = classes.map((I) => {
      let k = 1;
      let J = I;
      while (!J.isPrincipal()) {
        J = J.mul(I);
        k += 1;
      }
      return k;
    });
    const invariants = [];
    for (const { prime } of factorInteger(BigInt(h))) {
      const q = Number(prime);
      let before = 1;
      for (let k = 1, pk = q; before < h; k += 1, pk *= q) {
        const count = orders.filter((o) => pk % o === 0).length;
        if (count === before) break;
        const factors = Math.round(Math.log(count / before) / Math.log(q));
        // `factors` cyclic factors have order at least q^k
        for (let i = 0; i < factors; i += 1) {
          if (i < invariants.filter((v) => v.q === q).length) invariants.filter((v) => v.q === q)[i].e = k;
          else invariants.push({ q, e: k });
        }
        before = count;
      }
    }
    const cyc = invariants.map((v) => v.q ** v.e).sort((x, y) => x - y);
    return (this._classes = {
      classNumber: h,
      representatives: classes,
      orders,
      invariants: cyc,
      cyclic: new Set(invariants.map((v) => v.q)).size === invariants.length,
      text: cyc.length ? cyc.map((c) => `C${c}`).join(" × ") : "trivial",
      minkowskiBound: this.minkowskiBound,
      primesChecked: primes,
    });
  }
  classNumber() {
    return this.classGroup().classNumber;
  }
  /** every ideal is principal */
  get isPID() {
    return this.classNumber() === 1;
  }
  /** elements factor uniquely into irreducibles: for these rings, the same as being a PID */
  get isUFD() {
    return this.isPID;
  }
  /**
   * The factorisation of x into primes of the ring, when the ring has
   * unique factorisation: { unit, factors: [{ prime, exponent }] }.
   */
  factor(x) {
    x = this.of(x);
    if (x.isZero()) throw new RangeError("zero has no factorisation");
    if (!this.isUFD) throw new RangeError(`${this.name} does not have unique factorisation (its class number is ${this.classNumber()}); use factorizations(x) for the elements or ideal(x).factor() for ideals`);
    const out = [];
    let rest = x;
    for (const { prime: p } of factorInteger(abs(x.norm()))) {
      for (const P of this.primesAbove(p)) {
        const pi = this.canonical(P.generator());
        let k = 0;
        for (let q = rest.divExact(pi); q; q = rest.divExact(pi)) {
          rest = q;
          k += 1;
        }
        if (k) out.push({ prime: pi, exponent: k });
      }
    }
    return { unit: rest, factors: out };
  }
  toString() {
    return this.name;
  }
}

/** The ring of integers of Q(√d): quadraticRing(-5) is Z[√−5], quadraticRing(5) is Z[(1 + √5)/2]. */
export const quadraticRing = (d) => new QuadraticRing(d);
/** The Gaussian integers Z[i]. */
export const gaussianIntegers = () => new QuadraticRing(-1);
/** The Eisenstein integers Z[ω], ω = (1 + √−3)/2 (a sixth root of unity; ω − 1 is the usual cube root of unity). */
export const eisensteinIntegers = () => new QuadraticRing(-3);
