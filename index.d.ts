// Type declarations for @ratmath/structures.
import type { Rational } from "@ratmath/core";

/** The field interface shared across RatMath. Elements are whatever the field uses. */
export interface Field<E = any> {
  zero: E;
  one: E;
  add(a: E, b: E): E;
  sub(a: E, b: E): E;
  mul(a: E, b: E): E;
  div(a: E, b: E): E;
  neg(a: E): E;
  /** 0 when a is zero, otherwise nonzero. */
  sign(a: E): number | null;
  inv?(a: E): E;
  equals?(a: E, b: E): boolean;
  from?(x: unknown): E;
  of?(x: unknown): E;
  format?(a: E): string;
  name?: string;
}

export type Exponents = number[] | Record<string, number>;
export type MonomialOrder = "lex" | "grlex" | "grevlex" | ((a: number[], b: number[]) => number);
/** Anything a ring reads as one of its polynomials. C is the public coefficient type (Rational over Q). */
export type PolyLike<C = Rational> = MPoly<C> | string | number | bigint | Rational | C;

export interface RingOptions {
  field?: Field;
  order?: MonomialOrder;
}

export class PolyRing<C = Rational> {
  constructor(vars: string | string[], options?: RingOptions);
  readonly vars: readonly string[];
  readonly n: number;
  readonly field: Field | undefined;
  readonly order: "lex" | "grlex" | "grevlex" | "custom";
  readonly cmp: (a: number[], b: number[]) => number;
  readonly zero: MPoly<C>;
  readonly one: MPoly<C>;
  readonly gens: readonly MPoly<C>[];
  indexOf(v: string | number): number;
  gen(v: string | number): MPoly<C>;
  constant(c: unknown): MPoly<C>;
  monomial(exps: Exponents, c?: unknown): MPoly<C>;
  from(x: PolyLike<C>): MPoly<C>;
  poly(text: string): MPoly<C>;
  fromTerms(terms: Array<[unknown, Exponents]>): MPoly<C>;
  fromUnivariate(coeffs: unknown[], v?: string | number): MPoly<C>;
  withOrder(order: MonomialOrder): PolyRing<C>;
  equals(other: PolyRing<any>): boolean;
  toString(): string;
}

export function ring<C = Rational>(vars: string | string[], options?: RingOptions): PolyRing<C>;

export interface Term<C = Rational> {
  coefficient: C;
  exponents: number[];
  monomial: MPoly<C>;
}

export class MPoly<C = Rational> {
  readonly ring: PolyRing<C>;
  readonly length: number;
  terms(): Term<C>[];
  isZero(): boolean;
  isConstant(): boolean;
  isMonomial(): boolean;
  isHomogeneous(): boolean;
  constantTerm(): C;
  toConstant(): C;
  degree(v?: string | number): number;
  degrees(): Record<string, number>;
  variables(): string[];
  toUnivariate(): C[];
  leadingTerm(): MPoly<C>;
  leadingCoefficient(): C;
  leadingMonomial(): MPoly<C>;
  leadingExponents(): number[];
  coefficient(exps: Exponents): C;
  coefficients(v: string | number): MPoly<C>[];
  add(other: PolyLike<C>): MPoly<C>;
  sub(other: PolyLike<C>): MPoly<C>;
  mul(other: PolyLike<C>): MPoly<C>;
  neg(): MPoly<C>;
  scale(c: unknown): MPoly<C>;
  pow(n: number): MPoly<C>;
  equals(other: PolyLike<C>): boolean;
  monic(): MPoly<C>;
  divmod(divisor: PolyLike<C>): { quotient: MPoly<C>; remainder: MPoly<C> };
  divmod(divisors: PolyLike<C>[]): { quotients: MPoly<C>[]; remainder: MPoly<C> };
  rem(divisors: PolyLike<C> | PolyLike<C>[]): MPoly<C>;
  isDivisibleBy(g: PolyLike<C>): boolean;
  div(g: PolyLike<C>): MPoly<C>;
  diff(v: string | number, times?: number): MPoly<C>;
  gradient(): MPoly<C>[];
  subs(values: Record<string, PolyLike<C>>): MPoly<C>;
  at(point: unknown[] | Record<string, unknown>): C;
  compose<D = C>(values: Array<PolyLike<D>> | Record<string, PolyLike<D>>, target?: PolyRing<D>): MPoly<D>;
  homogenize(name?: string): MPoly<C>;
  into<D>(target: PolyRing<D>): MPoly<D>;
  toString(): string;
  toText(): string;
  toTeX(): string;
  toJSON(): { vars: string[]; terms: Array<[string, number[]]> };
}

/** The rationals, as the field object the library uses by default. */
export const QQ: Field<[bigint, bigint]> & { name: "Q" };
export const ORDERS: Readonly<Record<"lex" | "grlex" | "grevlex", (a: number[], b: number[]) => number>>;

export interface Factorisation<C = Rational, U = Rational> {
  unit: U;
  factors: Array<{ factor: MPoly<C>; multiplicity: number }>;
}

export function gcd<C>(first: PolyLike<C> | PolyLike<C>[], ...rest: PolyLike<C>[]): MPoly<C>;
export function lcm<C>(first: PolyLike<C> | PolyLike<C>[], ...rest: PolyLike<C>[]): MPoly<C>;
export function content<C>(f: MPoly<C>, v: string | number): MPoly<C>;
export function primitivePart<C>(f: MPoly<C>, v: string | number): MPoly<C>;
export function squareFree<C>(f: MPoly<C>): Factorisation<C, C>;
export function squareFreePart<C>(f: MPoly<C>): MPoly<C>;

export function det<C>(matrix: PolyLike<C>[][], ring?: PolyRing<C>): MPoly<C>;
export function sylvester<C>(f: MPoly<C>, g: PolyLike<C>, v: string | number): MPoly<C>[][];
export function resultant<C>(f: MPoly<C>, g: PolyLike<C>, v: string | number): MPoly<C>;
export function discriminant<C>(f: MPoly<C>, v: string | number): MPoly<C>;
export function jacobian<C>(list: MPoly<C>[]): MPoly<C>[][];
export function hessian<C>(f: MPoly<C>): MPoly<C>[][];

export function factor(f: MPoly<Rational>): Factorisation<Rational, Rational>;
export function isIrreducible(f: MPoly<Rational>): boolean;
export function integerPrimitive(f: MPoly<Rational>): { content: Rational; primitive: MPoly<Rational> };

export type FractionLike<C = Rational> = RationalFunction<C> | PolyLike<C>;

export class RationalFunction<C = Rational> {
  readonly field: FractionField<C>;
  readonly ring: PolyRing<C>;
  readonly num: MPoly<C>;
  readonly den: MPoly<C>;
  isZero(): boolean;
  isPolynomial(): boolean;
  toPolynomial(): MPoly<C>;
  add(other: FractionLike<C>): RationalFunction<C>;
  sub(other: FractionLike<C>): RationalFunction<C>;
  mul(other: FractionLike<C>): RationalFunction<C>;
  div(other: FractionLike<C>): RationalFunction<C>;
  neg(): RationalFunction<C>;
  inv(): RationalFunction<C>;
  pow(n: number): RationalFunction<C>;
  equals(other: FractionLike<C>): boolean;
  diff(v: string | number): RationalFunction<C>;
  subs(values: Record<string, PolyLike<C>>): RationalFunction<C>;
  at(point: unknown[] | Record<string, unknown>): C;
  toString(): string;
  toText(): string;
  toTeX(): string;
}

/** The field of fractions of a polynomial ring, as a RatMath field object. */
export interface FractionField<C = Rational> extends Field<RationalFunction<C>> {
  readonly ring: PolyRing<C>;
  readonly name: string;
  readonly gens: readonly RationalFunction<C>[];
  of(num: PolyLike<C>, den?: PolyLike<C>): RationalFunction<C>;
  from(x: FractionLike<C>): RationalFunction<C>;
  parse(text: string): RationalFunction<C>;
  inv(a: FractionLike<C>): RationalFunction<C>;
  pow(a: FractionLike<C>, n: number): RationalFunction<C>;
  equals(a: FractionLike<C>, b: FractionLike<C>): boolean;
  format(a: FractionLike<C>): string;
}

export function fractions<C = Rational>(R: PolyRing<C>): FractionField<C>;

export function elementary<C>(R: PolyRing<C>, k: number): MPoly<C>;
export function powerSum<C>(R: PolyRing<C>, k: number): MPoly<C>;
export function completeHomogeneous<C>(R: PolyRing<C>, k: number): MPoly<C>;
export function isSymmetric(f: MPoly<any>): boolean;
export function toElementary<C>(f: MPoly<C>, options?: { prefix?: string }): { polynomial: MPoly<C>; ring: PolyRing<C> };

export { primeField, galoisField, NumberField } from "@ratmath/algebraic";

// ---- finite groups -------------------------------------------------------------------------------

import type { PermGroup } from "@ratmath/algebraic";
export { PermGroup, cycleString, fromCycles, cycleType } from "@ratmath/algebraic";

/** A permutation of {0, …, n−1}: p[i] is the image of i. */
export type Perm = number[];
/** A square matrix over F_p, entries 0 … p−1. */
export type Matrix = number[][];
/** An element given by value, by printed label ("r²s", "(1 2 3)", "-i"), or by element number as { index }. */
export type ElementLike<E> = E | string | { index: number };

export interface ConjugacyClass<E> {
  readonly representative: E;
  readonly elements: E[];
  readonly indices: number[];
  readonly size: number;
  /** The order of each element in the class. */
  readonly order: number;
}

export interface SubgroupClassInfo<E> {
  readonly group: Group<E>;
  readonly order: number;
  readonly index: number;
  readonly conjugates: number;
  readonly normal: boolean;
  readonly members: Group<E>[];
}

export interface CosetInfo<E> {
  readonly representative: E;
  readonly elements: E[];
  readonly indices: number[];
}

export interface CayleyTable<E> {
  readonly labels: string[];
  readonly elements: E[];
  /** table[i][j] is the position (in `elements`) of the product of elements i and j. */
  readonly table: number[][];
  toString(): string;
}

export interface GroupSpec<E> {
  generators: E[];
  mul: (a: E, b: E) => E;
  identity?: E;
  key?: (x: E) => string;
  format?: (x: E) => string;
  name?: string | null;
  maxOrder?: number;
  labels?: Record<string, E> | null;
}

export class Group<E = any> {
  static generate<E>(spec: GroupSpec<E>): Group<E>;
  static fromElements<E>(elements: E[], mul: (a: E, b: E) => E, options?: { key?: (x: E) => string; format?: (x: E) => string; name?: string | null; generators?: E[] | null; identity?: E; labels?: Record<string, E> | null }): Group<E>;
  static fromTable<L = number>(table: number[][], options?: { labels?: L[] | null; name?: string | null; check?: boolean }): Group<L>;

  readonly elements: readonly E[];
  readonly order: number;
  name: string | null;
  /** For a subgroup: the group it was taken in. */
  readonly parent: Group<E> | null;
  /** Element numbers of the generators. */
  readonly gens: Int32Array;
  readonly generators: E[];
  readonly identity: E;
  /** Set on the named groups: their natural action (vertices, points, vectors). */
  natural?: Action<E, any>;
  /** Permutation groups: the degree, and the same group as algebraic's PermGroup. */
  degree?: number;
  toPermGroup?: () => PermGroup;
  /** Direct products: the factors. */
  factors?: Group<any>[];
  /** Quotients G/N: the coset of an element of G, the group G and the subgroup N. */
  cosetOf?: (x: any) => E;
  over?: Group<any>;
  normal?: Group<any>;

  indexOf(x: ElementLike<E>): number;
  element(x: ElementLike<E>): E;
  has(x: ElementLike<E>): boolean;
  key(x: E): string;
  format(x: E): string;
  mul(...xs: ElementLike<E>[]): E;
  inv(x: ElementLike<E>): E;
  pow(x: ElementLike<E>, k: number | bigint): E;
  /** g·x·g⁻¹. */
  conjugate(x: ElementLike<E>, g: ElementLike<E>): E;
  /** [a, b] = a⁻¹b⁻¹ab. */
  commutator(a: ElementLike<E>, b: ElementLike<E>): E;
  equal(a: ElementLike<E>, b: ElementLike<E>): boolean;
  mulIndex(i: number, j: number): number;
  invIndex(i: number): number;
  powIndex(i: number, k: number): number;
  orderIndex(i: number): number;

  elementOrder(x: ElementLike<E>): number;
  orderStatistics(): Map<number, number>;
  readonly exponent: number;
  readonly isAbelian: boolean;
  readonly isCyclic: boolean;
  readonly isTrivial: boolean;
  readonly isPGroup: boolean;
  readonly isSolvable: boolean;
  readonly isNilpotent: boolean;
  readonly nilpotencyClass: number | null;
  readonly isSimple: boolean;
  readonly isPerfect: boolean;
  cayleyTable(options?: { order?: ElementLike<E>[] | null }): CayleyTable<E>;

  subgroup(generators: ElementLike<E>[], options?: { name?: string | null }): Group<E>;
  trivialSubgroup(): Group<E>;
  subgroupFromIndices(list: number[], gens?: number[] | null, name?: string | null): Group<E>;
  closeIndices(gens: number[]): { list: number[]; inside: Uint8Array };
  maskOf(H: Group<E>): Uint8Array;
  parentIndices(): number[];
  toParentIndex(i: number): number;
  isSubgroup(H: Group<E>): boolean;
  equals(H: Group<E>): boolean;
  contains(H: Group<E>): boolean;
  index(H: Group<E> | ElementLike<E>[]): number;
  intersection(H: Group<E> | ElementLike<E>[], K: Group<E> | ElementLike<E>[]): Group<E>;
  join(...Hs: Array<Group<E> | ElementLike<E>[]>): Group<E>;
  product(H: Group<E> | ElementLike<E>[], K: Group<E> | ElementLike<E>[]): E[];

  centre(): Group<E>;
  center(): Group<E>;
  centraliser(S: ElementLike<E> | ElementLike<E>[] | Group<E>): Group<E>;
  centralizer(S: ElementLike<E> | ElementLike<E>[] | Group<E>): Group<E>;
  normaliser(H: Group<E> | ElementLike<E>[]): Group<E>;
  normalizer(H: Group<E> | ElementLike<E>[]): Group<E>;
  isNormal(H: Group<E> | ElementLike<E>[]): boolean;
  conjugateSubgroup(H: Group<E> | ElementLike<E>[], g: ElementLike<E>): Group<E>;
  normalClosure(S: ElementLike<E>[] | Group<E>): Group<E>;

  conjugacyClasses(): readonly ConjugacyClass<E>[];
  classEquation(): { order: number; centre: number; terms: number[]; toString(): string };

  derivedSubgroup(): Group<E>;
  commutatorSubgroup(): Group<E>;
  commutatorOf(H: Group<E> | ElementLike<E>[], K: Group<E> | ElementLike<E>[]): Group<E>;
  derivedSeries(): Group<E>[];
  lowerCentralSeries(): Group<E>[];
  upperCentralSeries(): Group<E>[];
  compositionSeries(): { series: Group<E>[]; factors: number[]; toString(): string };

  subgroups(options?: { maxOrder?: number }): readonly Group<E>[];
  subgroupClasses(options?: { maxOrder?: number }): readonly SubgroupClassInfo<E>[];
  normalSubgroups(): readonly Group<E>[];

  cosets(H: Group<E> | ElementLike<E>[], options?: { side?: "left" | "right" }): CosetInfo<E>[];
  quotient(N: Group<E> | ElementLike<E>[], options?: { name?: string | null }): Group<Coset>;

  sylowSubgroup(p: number): Group<E>;
  sylowSubgroups(p: number): Group<E>[];
  sylowCount(p: number): number;

  spanningTree(): { order: number[]; from: Int32Array; by: Int32Array };
  describe(): string;
  toString(): string;
}

/** An element of a quotient group: a coset gN. */
export class Coset {
  readonly id: number;
  /** Element number of the representative, in the group the quotient was taken of. */
  readonly repIndex: number;
  readonly indices: number[];
  readonly text: string;
  toString(): string;
}

/** Dihedral elements rᵏsᶠ and dicyclic elements aᵏxᶠ are pairs [k, f]. */
export type PowerPair = [number, 0 | 1];

export function cyclic(n: number): Group<number>;
export function dihedral(n: number): Group<PowerPair>;
export function dicyclic(n: number): Group<PowerPair>;
export function quaternion(): Group<PowerPair>;
export function kleinFour(): Group<[number, number]>;
export function permutationGroup(generators: Array<Perm | string> | PermGroup, n?: number, options?: { name?: string | null; maxOrder?: number }): Group<Perm>;
export function fromPermGroup(P: PermGroup, options?: { name?: string | null; maxOrder?: number }): Group<Perm>;
export function symmetric(n: number): Group<Perm>;
export function alternating(n: number): Group<Perm>;
export function matrixGroup(generators: Matrix[], p: number, options?: { name?: string | null; maxOrder?: number }): Group<Matrix>;
export function GL(n: number, p: number): Group<Matrix>;
export function SL(n: number, p: number): Group<Matrix>;
export function PSL(n: number, p: number): Group<Coset>;
export function heisenberg(p: number): Group<Matrix>;
export function directProduct(...groups: Array<Group<any> | { name?: string }>): Group<any[]>;
export function semidirectProduct<N, H>(N: Group<N>, H: Group<H>, act: (h: H) => (n: N) => N, options?: { name?: string | null }): Group<[N, H]>;
export function metacyclic(m: number, k: number, r: number, options?: { name?: string | null }): Group<[number, number]>;

export class Homomorphism<A = any, B = any> {
  constructor(source: Group<A>, target: Group<B>, map: ArrayLike<number>, options?: { name?: string | null });
  readonly source: Group<A>;
  readonly target: Group<B>;
  name: string | null;
  apply(x: ElementLike<A>): B;
  applyIndex(i: number): number;
  readonly table: number[];
  kernel(): Group<A>;
  image(): Group<B>;
  preimage(K: Group<B>): Group<A>;
  readonly isInjective: boolean;
  readonly isSurjective: boolean;
  readonly isIsomorphism: boolean;
  /** First this map, then g. */
  then<C>(g: Homomorphism<B, C>): Homomorphism<A, C>;
  inverse(): Homomorphism<B, A>;
  toString(): string;
}

export type HomSpec<A, B> = ((x: A) => B) | Array<ElementLike<B>> | Map<ElementLike<A>, ElementLike<B>>;
export function homomorphism<A, B>(G: Group<A>, H: Group<B>, f: HomSpec<A, B>, options?: { name?: string | null }): Homomorphism<A, B>;
export function tryHomomorphism<A, B>(G: Group<A>, H: Group<B>, f: HomSpec<A, B>, options?: { name?: string | null }): Homomorphism<A, B> | null;
export function identityMap<A>(G: Group<A>): Homomorphism<A, A>;
export function trivialMap<A, B>(G: Group<A>, H: Group<B>): Homomorphism<A, B>;
export function inclusion<A>(H: Group<A>, G: Group<A>): Homomorphism<A, A>;
export function projection(Q: Group<Coset>): Homomorphism<any, Coset>;
export function isomorphism<A, B>(G: Group<A>, H: Group<B>): Homomorphism<A, B> | null;
export function isIsomorphic(G: Group<any>, H: Group<any>): boolean;
export function automorphisms<A>(G: Group<A>, options?: { max?: number }): Homomorphism<A, A>[];
export function automorphismGroup<A>(G: Group<A>, options?: { max?: number }): Group<Homomorphism<A, A>> & { inner: Group<Homomorphism<A, A>>; of: Group<A> };
export function innerAutomorphism<A>(G: Group<A>, g: ElementLike<A>): Homomorphism<A, A>;
export function firstIsomorphism<A, B>(f: Homomorphism<A, B>): { kernel: Group<A>; quotient: Group<Coset>; image: Group<B>; iso: Homomorphism<Coset, B> };
export function secondIsomorphism<E>(G: Group<E>, H: Group<E>, N: Group<E>): { HN: Group<E>; HcapN: Group<E>; left: Group<Coset>; right: Group<Coset>; iso: Homomorphism<Coset, Coset> };
export function thirdIsomorphism<E>(G: Group<E>, N: Group<E>, K: Group<E>): { GN: Group<Coset>; KN: Group<Coset>; left: Group<Coset>; right: Group<Coset>; iso: Homomorphism<Coset, Coset> };

export interface OrbitCount<E> {
  readonly count: number;
  /** Σ_g |Fix(g)|. */
  readonly total: bigint;
  readonly order: number;
  readonly byClass: { representative: E; size: number; fixed: number }[];
  toString(): string;
}

export class Action<E = any, X = any> {
  constructor(group: Group<E>, points: X[], act: (g: E, x: X) => X, options?: { key?: (x: X) => string; format?: ((x: X) => string) | null; name?: string | null; check?: boolean });
  readonly group: Group<E>;
  readonly points: readonly X[];
  readonly act: (g: E, x: X) => X;
  readonly name: string | null;
  readonly degree: number;
  formatPoint(x: X): string;
  pointIndex(x: X): number;
  permutation(g: ElementLike<E>): number[];
  permutationIndex(i: number): Int32Array;
  apply(g: ElementLike<E>, x: X): X;
  orbit(x: X): X[];
  orbits(): X[][];
  stabiliser(x: X): Group<E>;
  stabilizer(x: X): Group<E>;
  fixedPoints(g: ElementLike<E>): X[];
  kernel(): Group<E>;
  readonly isTransitive: boolean;
  readonly isFaithful: boolean;
  readonly isFree: boolean;
  countOrbits(): OrbitCount<E>;
  toPermGroup(): PermGroup;
  onSubsets(k: number): Action<E, X[]>;
  onColourings<C = number>(colours: number | C[]): Action<E, C[]>;
  onColorings<C = number>(colours: number | C[]): Action<E, C[]>;
  toString(): string;
}

export function action<E, X>(group: Group<E>, points: X[], act: (g: E, x: X) => X, options?: { key?: (x: X) => string; format?: ((x: X) => string) | null; name?: string | null; check?: boolean }): Action<E, X>;
export function regularAction<E>(G: Group<E>): Action<E, E>;
export function conjugationAction<E>(G: Group<E>): Action<E, E>;
export function cosetAction<E>(G: Group<E>, H: Group<E>): Action<E, number> & { cosets: CosetInfo<E>[] };
export function burnside<E>(A: Action<E, any>): OrbitCount<E>;

export const SMALL_ORDER_LIMIT: 20;
export function smallGroups(n: number): Group<any>[];
export function numberOfGroups(n: number): number;
export function identify<E>(G: Group<E>): { name: string; group: Group<any>; iso: Homomorphism<E, any> } | null;

// ---- Gröbner bases ------------------------------------------------------------

export interface GroebnerOptions {
  /** A monomial order to use in place of the ring's own. */
  order?: MonomialOrder;
  /** Track each basis element as a combination of the generators while computing. */
  certificates?: boolean;
  /** Give up after this many S-pairs (default 200000). */
  maxPairs?: number;
}

export interface GroebnerStats {
  pairs: number;
  zero: number;
  skipped: number;
  added: number;
}

export class GroebnerBasis<C = Rational> {
  readonly ring: PolyRing<C>;
  readonly generators: MPoly<C>[];
  /** The reduced Gröbner basis: monic, sorted by leading monomial. */
  readonly basis: MPoly<C>[];
  readonly stats: GroebnerStats;
  readonly length: number;
  readonly leadingExponents: number[][];
  readonly leadingMonomials: MPoly<C>[];
  /** The ideal is the whole ring (the basis is {1}). */
  readonly isTrivial: boolean;
  readonly isZero: boolean;
  /** The remainder on division by the basis: the normal form of f modulo the ideal. */
  reduce(f: PolyLike<C>): MPoly<C>;
  divide(f: PolyLike<C>): { quotients: MPoly<C>[]; remainder: MPoly<C> };
  contains(f: PolyLike<C>): boolean;
  /** Each basis element as a combination of the generators (rows match basis). */
  readonly cofactors: MPoly<C>[][];
  /** f = Σ cofactors[i]·generators[i], for f in the ideal; throws otherwise. */
  express(f: PolyLike<C>): { cofactors: MPoly<C>[]; check: boolean };
  /** Buchberger's criterion, re-checked from scratch. */
  verify(): { generatorsReduce: boolean; sPairsReduce: boolean; basisInIdeal: boolean; sPairsChecked: number; ok: boolean };
  /** The dimension of the solution set; −1 for the whole ring. */
  readonly dimension: number;
  readonly isZeroDimensional: boolean;
  /** The monomials no leading monomial divides: a basis of the quotient (zero-dimensional ideals). */
  standardMonomials(options?: { limit?: number }): MPoly<C>[];
  /** How many standard monomials there are; Infinity when the ideal is not zero-dimensional. */
  readonly quotientDimension: number;
  toString(): string;
}

export function groebner<C = Rational>(polys: MPoly<C>[], options?: GroebnerOptions): GroebnerBasis<C>;
export function sPolynomial<C = Rational>(f: MPoly<C>, g: MPoly<C>): MPoly<C>;
/** A block order that puts the listed variable indices above all others, so a Gröbner basis eliminates them. */
export function eliminationOrder(n: number, eliminate: number[]): (a: number[], b: number[]) => number;

// ---- ideals and quotient rings ---------------------------------------------------

export type RadicalCertificate<C = Rational> =
  | { holds: false }
  | {
      holds: true;
      /** f^power = Σ cofactors[i]·generators[i] */
      power: number;
      cofactors: MPoly<C>[];
      check: boolean;
      extended: { ring: PolyRing<C>; generators: MPoly<C>[]; cofactors: MPoly<C>[] };
    };

export class Ideal<C = Rational> {
  constructor(ring: PolyRing<C>, generators: MPoly<C>[]);
  readonly ring: PolyRing<C>;
  readonly generators: MPoly<C>[];
  readonly groebner: GroebnerBasis<C>;
  readonly basis: MPoly<C>[];
  readonly isZero: boolean;
  readonly isWhole: boolean;
  readonly isProper: boolean;
  contains(f: PolyLike<C>): boolean;
  reduce(f: PolyLike<C>): MPoly<C>;
  express(f: PolyLike<C>): { cofactors: MPoly<C>[]; check: boolean };
  isSubsetOf(J: Ideal<C>): boolean;
  equals(J: Ideal<C>): boolean;
  add(J: Ideal<C> | PolyLike<C>): Ideal<C>;
  mul(J: Ideal<C> | PolyLike<C>): Ideal<C>;
  pow(n: number): Ideal<C>;
  intersect(J: Ideal<C> | PolyLike<C>): Ideal<C>;
  /** I : J, the polynomials that multiply J into I. */
  quotient(J: Ideal<C> | PolyLike<C>): Ideal<C>;
  /** I : f^∞ (or I : J^∞). */
  saturate(f: Ideal<C> | PolyLike<C>): Ideal<C>;
  /** Whether some power of f lies in the ideal (f vanishes wherever the generators do). */
  radicalContains(f: PolyLike<C>): boolean;
  radicalContains(f: PolyLike<C>, options: { certificate: true }): RadicalCertificate<C>;
  exponentOf(f: PolyLike<C>, options?: { limit?: number }): number | null;
  eliminate(vars: string | number | Array<string | number>): Ideal<C>;
  readonly dimension: number;
  readonly isZeroDimensional: boolean;
  /** dim of R/I as a vector space: solutions counted with multiplicity. */
  readonly degree: number;
  standardMonomials(options?: { limit?: number }): MPoly<C>[];
  /** Minimal polynomial of f modulo the ideal, constant first (zero-dimensional ideals). */
  minimalPolynomial(f: PolyLike<C>): C[];
  radical(): Ideal<C>;
  isRadical(): boolean;
  isMaximal(): boolean;
  /** null when the test implemented here cannot decide. */
  isPrime(): boolean | null;
  quotientRing(): QuotientRing<C>;
  solve(options?: SolveOptions): Solutions;
  toString(): string;
}

export function ideal<C = Rational>(ring: PolyRing<C>, generators: Array<PolyLike<C>>): Ideal<C>;
export function ideal<C = Rational>(generators: MPoly<C>[]): Ideal<C>;

export interface ChineseRemainder<C = Rational> {
  comaximal: boolean;
  /** a ∈ I and b ∈ J with a + b = 1 (when comaximal). */
  a?: MPoly<C>;
  b?: MPoly<C>;
  /** I + J, when the ideals are not comaximal. */
  sum?: Ideal<C>;
  intersection?: Ideal<C>;
  productIsIntersection?: boolean;
  /** A polynomial congruent to f mod I and to g mod J. */
  lift?(f: PolyLike<C>, g: PolyLike<C>): MPoly<C>;
}
export function chineseRemainder<C = Rational>(I: Ideal<C>, J: Ideal<C>): ChineseRemainder<C>;

export class QuotientElement<C = Rational> {
  readonly parent: QuotientRing<C>;
  /** The normal form: the representative no leading monomial divides. */
  readonly rep: MPoly<C>;
  add(o: QuotientElement<C> | PolyLike<C>): QuotientElement<C>;
  sub(o: QuotientElement<C> | PolyLike<C>): QuotientElement<C>;
  neg(): QuotientElement<C>;
  mul(o: QuotientElement<C> | PolyLike<C>): QuotientElement<C>;
  pow(n: number): QuotientElement<C>;
  equals(o: QuotientElement<C> | PolyLike<C>): boolean;
  isZero(): boolean;
  isUnit(): boolean;
  inverse(): QuotientElement<C>;
  div(o: QuotientElement<C> | PolyLike<C>): QuotientElement<C>;
  isZeroDivisor(): boolean;
  isNilpotent(): boolean;
  minimalPolynomial(): C[];
  /** Multiplication by this element on the standard monomials (finite quotients). */
  matrix(): C[][];
  toString(): string;
}

export class QuotientRing<C = Rational> {
  constructor(I: Ideal<C>);
  readonly ring: PolyRing<C>;
  readonly ideal: Ideal<C>;
  readonly gens: QuotientElement<C>[];
  readonly zero: QuotientElement<C>;
  readonly one: QuotientElement<C>;
  of(x: QuotientElement<C> | PolyLike<C>): QuotientElement<C>;
  /** Vector-space dimension over the field; Infinity when not finite. */
  readonly dimension: number;
  readonly isFinite: boolean;
  readonly basis: MPoly<C>[];
  readonly isZeroRing: boolean;
  isField(): boolean;
  isIntegralDomain(): boolean | null;
  isReduced(): boolean;
  /** The quotient as a Field, usable as the coefficients of a new ring (the ideal must be maximal). */
  asField(options?: { name?: string | null }): Field;
  toString(): string;
}

export interface SolveOptions {
  /** Also list the solutions that are not real. */
  complex?: boolean;
}

export interface Solution {
  real: boolean;
  /** The root of the component's minimal polynomial this solution belongs to. */
  theta: any;
  field: any;
  /** Each coordinate as an element of the number field Q(θ). */
  coordinates: Record<string, any>;
  /** Each coordinate as an exact algebraic number (computed on first use). */
  readonly values: Record<string, any>;
  /** Floating-point coordinates (real solutions). */
  approx?: Record<string, number>;
  toString(): string;
}

export interface SolutionComponent {
  field: any;
  minimalPolynomial: bigint[];
  degree: number;
  realCount: number;
  coordinates: Record<string, any>;
}

export interface Solutions {
  finite: boolean;
  /** 0 for finitely many solutions, −1 for none because the ideal is everything, otherwise the dimension of the solution set. */
  dimension: number;
  variables: string[];
  /** Distinct solutions over the complex numbers. */
  count?: number;
  /** Solutions counted with multiplicity. */
  multiplicity?: number;
  radical?: boolean;
  real?: Solution[];
  complex?: Solution[];
  components?: SolutionComponent[];
  parametrization?: {
    form: Record<string, number>;
    minimalPolynomial: Rational[];
    coordinates: Record<string, Rational[]>;
  };
}

export function solve(ring: PolyRing, polys: Array<PolyLike>, options?: SolveOptions): Solutions;
export function solve(polys: MPoly[], options?: SolveOptions): Solutions;
export function solve(ideal: Ideal, options?: SolveOptions): Solutions;

// ---- the integers mod n -------------------------------------------------------

export type IntLike = number | bigint;

/** Solve x ≡ residues[i] (mod moduli[i]); null when the congruences contradict each other. */
export function crt(residues: IntLike[], moduli: IntLike[]): { value: bigint; modulus: bigint } | null;

export interface ZnIdeal {
  generator: bigint;
  size: bigint;
  quotient: string;
  isMaximal: boolean;
  isPrime: boolean;
}

export class ZnRing {
  constructor(n: IntLike);
  readonly n: bigint;
  readonly characteristic: bigint;
  readonly size: bigint;
  of(a: IntLike): bigint;
  add(a: IntLike, b: IntLike): bigint;
  sub(a: IntLike, b: IntLike): bigint;
  mul(a: IntLike, b: IntLike): bigint;
  neg(a: IntLike): bigint;
  pow(a: IntLike, e: IntLike): bigint;
  isUnit(a: IntLike): boolean;
  inverse(a: IntLike): bigint;
  div(a: IntLike, b: IntLike): bigint;
  isZeroDivisor(a: IntLike): boolean;
  readonly isField: boolean;
  readonly isIntegralDomain: boolean;
  readonly isReduced: boolean;
  elements(): bigint[];
  units(): bigint[];
  readonly unitCount: bigint;
  zeroDivisors(): bigint[];
  nilpotents(): bigint[];
  idempotents(): bigint[];
  additiveOrder(a: IntLike): bigint;
  /** The multiplicative order of a unit. */
  orderOf(a: IntLike): bigint;
  unitGroup(): { order: bigint; invariants: bigint[]; cyclic: boolean; text: string };
  primitiveRoot(): bigint | null;
  ideals(): ZnIdeal[];
  /** Z/n as a product of rings Z/pᵏ, with the idempotent that picks out each factor. */
  decompose(): { moduli: bigint[]; idempotents: bigint[]; to(a: IntLike): bigint[]; from(residues: IntLike[]): bigint };
  table(op?: "add" | "mul"): bigint[][];
  /** The field Z/p (n must be prime), for use as polynomial coefficients. */
  asField(): Field;
  toString(): string;
}
export function Zn(n: IntLike): ZnRing;

// ---- quadratic integer rings ----------------------------------------------------

export type QuadLike = QuadraticInteger | IntLike;

export class QuadraticInteger {
  readonly ring: QuadraticRing;
  /** The element is a + bω. */
  readonly a: bigint;
  readonly b: bigint;
  add(o: QuadLike): QuadraticInteger;
  sub(o: QuadLike): QuadraticInteger;
  neg(): QuadraticInteger;
  mul(o: QuadLike): QuadraticInteger;
  /** Negative powers are allowed for units. */
  pow(k: IntLike): QuadraticInteger;
  conj(): QuadraticInteger;
  norm(): bigint;
  trace(): bigint;
  isZero(): boolean;
  equals(o: QuadLike): boolean;
  isUnit(): boolean;
  inverse(): QuadraticInteger;
  /** The quotient when o divides this exactly, otherwise null. */
  divExact(o: QuadLike): QuadraticInteger | null;
  isDivisibleBy(o: QuadLike): boolean;
  isAssociateOf(o: QuadLike): boolean;
  /** The sign of the real number (real quadratic rings). */
  sign(): number;
  toNumber(): number;
  toString(): string;
}

export class QuadraticIdeal {
  readonly ring: QuadraticRing;
  readonly norm: bigint;
  /** Two elements that generate the ideal. */
  readonly generators: QuadraticInteger[];
  readonly isWhole: boolean;
  contains(x: QuadLike): boolean;
  equals(o: QuadraticIdeal): boolean;
  isSubsetOf(o: QuadraticIdeal): boolean;
  add(o: QuadraticIdeal): QuadraticIdeal;
  mul(o: QuadraticIdeal): QuadraticIdeal;
  pow(k: number): QuadraticIdeal;
  conj(): QuadraticIdeal;
  /** A single generator, or null when the ideal is not principal. */
  generator(): QuadraticInteger | null;
  isPrincipal(): boolean;
  isPrime(): boolean;
  /** Same ideal class: the ideals differ by a principal factor. */
  isEquivalentTo(o: QuadraticIdeal): boolean;
  factor(): Array<{ prime: QuadraticIdeal; exponent: number }>;
  toString(): string;
}

export interface QuadraticUnits {
  finite: boolean;
  list?: QuadraticInteger[];
  count?: number;
  fundamental?: QuadraticInteger;
  norm?: bigint;
  text?: string;
}

export interface ClassGroup {
  classNumber: number;
  representatives: QuadraticIdeal[];
  /** The order of each representative's class. */
  orders: number[];
  /** Prime-power cyclic factors. */
  invariants: number[];
  cyclic: boolean;
  text: string;
  minkowskiBound: number;
  primesChecked: QuadraticIdeal[];
}

export class QuadraticRing {
  /** The ring of integers of Q(√d), for square-free d other than 0 and 1. */
  constructor(d: IntLike);
  readonly d: bigint;
  readonly discriminant: bigint;
  readonly zero: QuadraticInteger;
  readonly one: QuadraticInteger;
  /** ω = √d, or (1 + √d)/2 when d ≡ 1 mod 4. */
  readonly omega: QuadraticInteger;
  readonly isReal: boolean;
  /** a + bω. */
  of(a: QuadLike, b?: IntLike): QuadraticInteger;
  /** x + y√d (which must be an integer of the ring). */
  sqrt(x: IntLike, y: IntLike): QuadraticInteger;
  readonly units: QuadraticUnits;
  readonly fundamentalUnit: QuadraticInteger | null;
  readonly isNormEuclidean: boolean;
  divmod(a: QuadLike, b: QuadLike): { quotient: QuadraticInteger; remainder: QuadraticInteger; euclidean: boolean };
  gcd(a: QuadLike, b: QuadLike): { gcd: QuadraticInteger; s: QuadraticInteger; t: QuadraticInteger; steps: unknown[] };
  /** The preferred associate of x. */
  canonical(x: QuadLike): QuadraticInteger;
  /** Elements of norm m (−m too with bothSigns), one from each class of associates. */
  elementsOfNorm(m: IntLike, options?: { bothSigns?: boolean; limit?: number }): QuadraticInteger[];
  isIrreducible(x: QuadLike): boolean;
  isPrimeElement(x: QuadLike): boolean;
  irreducibleDivisors(x: QuadLike): QuadraticInteger[];
  /** Every way to write x as a unit times irreducibles, up to order and associates. */
  factorizations(x: QuadLike, options?: { limit?: number }): Array<{ unit: QuadraticInteger; factors: QuadraticInteger[] }>;
  ideal(...generators: QuadLike[]): QuadraticIdeal;
  splitting(p: IntLike): { prime: IntLike; kind: "split" | "inert" | "ramified"; primes: QuadraticIdeal[]; principal: boolean };
  primesAbove(p: IntLike): QuadraticIdeal[];
  readonly minkowskiBound: number;
  classGroup(): ClassGroup;
  classNumber(): number;
  readonly isPID: boolean;
  readonly isUFD: boolean;
  /** The factorisation into primes; throws when the ring does not have unique factorisation. */
  factor(x: QuadLike): { unit: QuadraticInteger; factors: Array<{ prime: QuadraticInteger; exponent: number }> };
  toString(): string;
}
export function quadraticRing(d: IntLike): QuadraticRing;
export function gaussianIntegers(): QuadraticRing;
export function eisensteinIntegers(): QuadraticRing;
