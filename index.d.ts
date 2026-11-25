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
