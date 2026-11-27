import { describe, expect, test } from "bun:test";
import { cosMinimal } from "@ratmath/algebraic";
import { sign, isEven, isOdd, inversions, fromCycles, symmetric, alternating, constructible, constructiblePolygon, regularPolygon } from "../index.js";

describe("the sign of a permutation", () => {
  test("cycles of each length", () => {
    expect(sign("(1 2 3 4)", 4)).toBe(-1);
    expect(isOdd("(1 2 3 4)", 4)).toBe(true);
    expect(sign("(1 2 3)", 3)).toBe(1);
    expect(sign("(1 2)(3 4)", 4)).toBe(1);
    expect(sign("()", 3)).toBe(1);
    expect(sign([1, 0, 2])).toBe(-1);
    expect(isEven(fromCycles("(1 2 3)(4 5 6 7)", 7))).toBe(false);
    expect(sign("(1 5)")).toBe(-1); // n read from the largest point
    expect(() => sign([0, 0, 1])).toThrow(RangeError);
  });

  test("the sign is the parity of the inversions (independent check, all of S5)", () => {
    const perms = (n) => (n === 0 ? [[]] : perms(n - 1).flatMap((p) => Array.from({ length: n }, (_, i) => [...p.slice(0, i), n - 1, ...p.slice(i)])));
    let even = 0;
    for (const p of perms(5)) {
      expect(sign(p)).toBe(inversions(p) % 2 === 0 ? 1 : -1);
      if (isEven(p)) even += 1;
    }
    expect(even).toBe(60);
  });

  test("agrees with membership in the alternating group", () => {
    const A4 = alternating(4);
    for (const text of ["(1 2 3 4)", "(1 2 3)", "(1 2)(3 4)", "(1 3)", "()"]) {
      expect(A4.has(fromCycles(text, 4))).toBe(isEven(text, 4));
    }
  });

  test("sign is a homomorphism", () => {
    const S = symmetric(4);
    const els = S.elements ?? [...Array(S.order).keys()].map((i) => S.element(i));
    for (const a of els.slice(0, 24)) for (const b of els.slice(0, 24)) expect(sign(S.mul(a, b))).toBe(sign(a) * sign(b));
  });
});

describe("Gauss–Wantzel in constructible", () => {
  test("the 17-gon: cos(2π/17) has degree 8 and is constructible", () => {
    const c = constructible(cosMinimal(17));
    expect(c.constructible).toBe(true);
    expect(c.factors[0].cosine).toEqual({ n: 17, scale: 1 });
    expect(c.reason).toContain("Gauss–Wantzel");
  });

  test("2cos(2π/17) is recognised too", () => {
    const m = cosMinimal(17).map((c, k) => `${c}/${2n ** BigInt(k)}`);
    const c = constructible(m);
    expect(c.constructible).toBe(true);
    expect(c.factors[0].cosine).toEqual({ n: 17, scale: 2 });
  });

  test("agrees with regularPolygon for every n with φ(n)/2 a power of 2 above 7", () => {
    for (const n of [17, 32, 34, 40, 48, 60, 51, 64, 68, 80]) {
      const c = constructible(cosMinimal(n));
      expect(c.constructible).toBe(regularPolygon(n).constructible);
      expect(c.constructible).toBe(true);
    }
  });

  test("other degree-8 polynomials stay undecided; non powers of 2 are refused", () => {
    expect(constructible("x^8 - 2").constructible).toBe(null);
    expect(constructible(cosMinimal(19)).constructible).toBe(false);
  });

  test("constructiblePolygon is regularPolygon", () => {
    for (const n of [3, 5, 7, 9, 15, 17, 51, 255, 257]) expect(constructiblePolygon(n).constructible).toBe(regularPolygon(n).constructible);
    expect(constructiblePolygon(7).constructible).toBe(false);
  });
});
