import { describe, expect, test } from "bun:test";
import { crt, Zn, quadraticRing, gaussianIntegers, eisensteinIntegers } from "../index.js";

const strs = (xs) => xs.map(String);

describe("the integers mod n", () => {
  const Z12 = Zn(12);
  test("units, zero divisors, nilpotents, idempotents", () => {
    expect(strs(Z12.units())).toEqual(["1", "5", "7", "11"]);
    expect(strs(Z12.zeroDivisors())).toEqual(["2", "3", "4", "6", "8", "9", "10"]);
    expect(strs(Z12.nilpotents())).toEqual(["0", "6"]);
    expect(strs(Z12.idempotents())).toEqual(["0", "1", "4", "9"]);
    expect(Z12.isField).toBe(false);
    expect(Zn(7).isField).toBe(true);
    expect(Z12.isReduced).toBe(false);
    expect(Zn(30).isReduced).toBe(true);
  });
  test("arithmetic and inverses", () => {
    const Z7 = Zn(7);
    expect(Z7.inverse(3)).toBe(5n);
    expect(Z7.div(1, 3)).toBe(5n);
    expect(Z7.pow(3, 6)).toBe(1n);
    expect(Z7.pow(3, -1)).toBe(5n);
    expect(() => Z12.inverse(4)).toThrow();
    expect(Z12.add(7, 8)).toBe(3n);
    expect(Z12.neg(5)).toBe(7n);
  });
  test("the unit group and primitive roots", () => {
    expect(Zn(7).primitiveRoot()).toBe(3n);
    expect(Zn(8).primitiveRoot()).toBe(null);
    expect(Z12.unitGroup().text).toBe("C2 × C2");
    expect(Zn(15).orderOf(2)).toBe(4n);
    expect(Zn(15).additiveOrder(6)).toBe(5n);
  });
  test("ideals and the decomposition into prime-power parts", () => {
    expect(Z12.ideals().length).toBe(6);
    expect(Z12.ideals().filter((I) => I.isMaximal).map((I) => I.generator)).toEqual([2n, 3n]);
    const d = Z12.decompose();
    expect(d.moduli).toEqual([4n, 3n]);
    expect(d.idempotents).toEqual([9n, 4n]);
  });
  test("the Chinese remainder theorem", () => {
    expect(crt([2, 3, 2], [3, 5, 7])).toEqual({ value: 23n, modulus: 105n });
    expect(crt([1, 3], [4, 6])).toEqual({ value: 9n, modulus: 12n });
    expect(crt([1, 2], [4, 6])).toBe(null);
  });
});

describe("quadratic integer rings", () => {
  test("class numbers", () => {
    const h = (d) => quadraticRing(d).classNumber();
    expect([-1, -2, -3, -7, -11, -19, -43, -67, -163].map(h)).toEqual([1, 1, 1, 1, 1, 1, 1, 1, 1]);
    expect([-5, -6, -10, -13, -15].map(h)).toEqual([2, 2, 2, 2, 2]);
    expect([-23, -14, -47, 10, 79, 82, 229, 401].map(h)).toEqual([3, 4, 5, 2, 3, 4, 3, 5]);
    expect(quadraticRing(-65).classGroup().text).toBe("C2 × C4");
    expect(quadraticRing(-21).classGroup().text).toBe("C2 × C2");
    expect(quadraticRing(-14).classGroup().cyclic).toBe(true);
    expect(() => quadraticRing(12)).toThrow("square-free");
  });
  test("fundamental units", () => {
    const u = (d) => String(quadraticRing(d).fundamentalUnit);
    expect(u(2)).toBe("1 + √2");
    expect(u(3)).toBe("2 + √3");
    expect(u(5)).toBe("(1 + √5)/2");
    expect(u(7)).toBe("8 + 3√7");
    expect(u(13)).toBe("(3 + √13)/2");
    expect(u(94)).toBe("2143295 + 221064√94");
    expect(quadraticRing(2).fundamentalUnit.norm()).toBe(-1n);
    expect(quadraticRing(-1).units.count).toBe(4);
    expect(eisensteinIntegers().units.count).toBe(6);
    expect(quadraticRing(-5).units.count).toBe(2);
  });
  test("Gaussian integers: division, gcd, factoring", () => {
    const G = gaussianIntegers();
    const a = G.of(11, 3);
    const b = G.of(1, 8);
    const { quotient, remainder } = G.divmod(a, b);
    expect(quotient.mul(b).add(remainder).equals(a)).toBe(true);
    expect(remainder.norm() < b.norm()).toBe(true);
    const g = G.gcd(a, b);
    expect(String(g.gcd)).toBe("2 + i");
    expect(g.s.mul(a).add(g.t.mul(b)).equals(g.gcd)).toBe(true);
    const f = G.factor(G.of(30));
    expect(f.factors.map(({ prime, exponent }) => `${prime}^${exponent}`)).toEqual(["1 + i^2", "3^1", "2 − i^1", "2 + i^1"]);
    const back = f.factors.reduce((p, { prime, exponent }) => p.mul(prime.pow(exponent)), f.unit);
    expect(back.equals(G.of(30))).toBe(true);
    expect(G.isIrreducible(G.of(3))).toBe(true);
    expect(G.isIrreducible(G.of(5))).toBe(false);
    expect(G.isUFD).toBe(true);
    expect(G.isNormEuclidean).toBe(true);
  });
  test("Z[√−5]: 6 factors in two ways, and ideals repair it", () => {
    const K = quadraticRing(-5);
    expect(K.isUFD).toBe(false);
    const fs = K.factorizations(K.of(6)).map((f) => strs(f.factors));
    expect(fs).toEqual([["2", "3"], ["1 + √−5", "1 − √−5"]]);
    expect(K.isIrreducible(K.of(2))).toBe(true);
    expect(K.isPrimeElement(K.of(2))).toBe(false);
    expect(() => K.factor(K.of(6))).toThrow("unique factorisation");
    const parts = K.ideal(K.of(6)).factor().map(({ prime, exponent }) => `${prime}^${exponent}`);
    expect(parts).toEqual(["(2, 1 + √−5)^2", "(3, 2 + √−5)^1", "(3, 1 + √−5)^1"]);
    const P = K.ideal(K.of(2), K.sqrt(1, 1));
    expect(P.isPrincipal()).toBe(false);
    expect(P.isPrime()).toBe(true);
    expect(P.mul(P).equals(K.ideal(K.of(2)))).toBe(true);
    expect(String(P.mul(P).generator())).toBe("2");
  });
  test("how primes split", () => {
    const K = quadraticRing(-5);
    expect([2, 3, 5, 7, 11].map((p) => K.splitting(p).kind)).toEqual(["ramified", "split", "ramified", "split", "inert"]);
    const G = gaussianIntegers();
    expect([2, 3, 5, 7, 13].map((p) => G.splitting(p).kind)).toEqual(["ramified", "inert", "split", "inert", "split"]);
  });
  test("elements of a given norm", () => {
    expect(strs(gaussianIntegers().elementsOfNorm(5))).toEqual(["2 + i", "2 − i"]);
    expect(quadraticRing(-5).elementsOfNorm(2)).toEqual([]);
    expect(strs(quadraticRing(2).elementsOfNorm(7))).toEqual(["3 + √2", "5 + 3√2"]);
    expect(strs(quadraticRing(2).elementsOfNorm(-1))).toEqual(["1 + √2"]);
    expect(quadraticRing(3).elementsOfNorm(-1)).toEqual([]);
  });
  test("arithmetic of elements", () => {
    const E = eisensteinIntegers();
    const w = E.omega;
    expect(w.pow(6).equals(E.one)).toBe(true);
    expect(w.pow(3).equals(E.one.neg())).toBe(true);
    expect(w.norm()).toBe(1n);
    expect(w.inverse().mul(w).equals(E.one)).toBe(true);
    const K = quadraticRing(2);
    const e = K.fundamentalUnit;
    expect(e.pow(-2).mul(e.pow(2)).equals(K.one)).toBe(true);
    expect(K.of(3, 1).divExact(K.of(1, 1)) === null).toBe(false);
    expect(K.of(1, 1).isAssociateOf(K.of(3, 2))).toBe(true);
  });
});
