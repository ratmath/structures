import { describe, expect, test } from "bun:test";
import { solveByRadicals, rootOfUnityByRadicals, regularPolygon, constructible, classicalProblems } from "../index.js";

const first = (p, k = 0) => solveByRadicals(p).factors[0].roots[k].expression;
const allAgree = (p) => solveByRadicals(p).factors.every((f) => f.roots.every((r) => r.expression.check(30).agrees === true));

describe("solving by radicals", () => {
  test("quadratics give the quadratic formula", () => {
    expect(String(first("x^2 - x - 1", 1))).toBe("(1 + √5)/2");
    expect(String(first("x^2 + x + 1"))).toBe("(−1 + i·√3)/2");
    expect(String(first("x^2 - 8", 1))).toBe("2·√2");
    expect(String(first("3x^2 - 2", 1))).toBe("√6/3");
  });
  test("cubics give Cardano's formula", () => {
    expect(String(first("x^3 - 2"))).toBe("∛2");
    expect(String(first("x^3 + x + 1"))).toBe("(−∛(108 + 12·√93) − ∛(108 − 12·√93))/6");
    // three real roots: the radicals have to pass through complex numbers
    const e = first("x^3 - 3x + 1");
    expect(String(e)).toBe("(−∛(4 + 4·i·√3) − ∛(4 − 4·i·√3))/2");
    expect(e.evaluate(20).real).toBe(true);
    expect(e.approx()).toBeCloseTo(-1.8793852415718169, 12);
    expect(e.indices).toEqual([2, 3]);
  });
  test("every root of a solvable polynomial is written and checked", () => {
    for (const p of ["x^3 - 2", "x^4 - 2", "x^4 + x + 1", "x^4 - 10x^2 + 1", "x^4 + 1", "2x^3 - 4x + 1", "x^6 - 3x^2 - 1", "x^6 + x^3 + 1"]) expect(allAgree(p)).toBe(true);
  });
  test("solvable quintics, by Lagrange resolvents", () => {
    const r = solveByRadicals("x^5 - 5x + 12");
    expect(r.solvable).toBe(true);
    expect(r.factors[0].group).toBe("D5");
    expect(r.factors[0].chain).toEqual([2, 5]);
    expect(r.factors[0].roots.every((x) => x.expression.check(40).agrees)).toBe(true);
    expect(r.factors[0].roots[0].expression.indices).toEqual([2, 5]);
    expect(String(first("x^5 - 2"))).toBe("⁵√2");
    expect(allAgree("x^5 + x^4 - 4x^3 - 3x^2 + 3x + 1")).toBe(true);
  });
  test("degree 7", () => {
    const r = solveByRadicals("x^7 + 7x^3 + 7x^2 + 7x - 1");
    expect(r.factors[0].group).toBe("D7");
    expect(r.factors[0].roots[0].expression.check(30).agrees).toBe(true);
    expect(r.factors[0].roots[0].expression.definitions.some((d) => d.name === "ζ₇")).toBe(true);
  });
  test("an unsolvable quintic is reported with its group", () => {
    const r = solveByRadicals("x^5 - x - 1");
    expect(r.solvable).toBe(false);
    expect(r.factors[0].group).toBe("S5");
    expect(r.factors[0].roots).toBe(null);
    expect(r.reason).toMatch(/not solvable/);
    expect(solveByRadicals("x^5 + 20x + 16").factors[0].group).toBe("A5");
  });
  test("reducible polynomials are solved factor by factor", () => {
    const r = solveByRadicals("(x^2 - 2)(x^3 - 5)(x - 3)".replace(/\)\(/g, ")*("));
    expect(r.factors.map((f) => f.degree)).toEqual([1, 2, 3]);
    expect(String(r.factors[0].roots[0].expression)).toBe("3");
    expect(String(r.factors[0].roots[0].expression.rational)).toBe("3");
    expect(r.solvable).toBe(true);
  });
  test("definitions, values and TeX", () => {
    const e = first("x^4 + x + 1");
    expect(e.inline({ maxLength: 60 })).toBe(null);
    expect(e.definitions.length).toBeGreaterThan(3);
    expect(e.definitions[0].kind).toBe("root");
    expect(e.definitions[0].index).toBe(3);
    const v = e.evaluate(30);
    expect(v.re.startsWith("−0.72713608449119683997667565")).toBe(true);
    expect(first("x^2 - x - 1", 1).toTeX()).toBe("\\frac{1 + \\sqrt{5}}{2}");
    expect(() => solveByRadicals("x^8 - 2")).toThrow(/degree 8/);
    expect(() => solveByRadicals("x^2 - y")).toThrow();
  });
});

describe("roots of unity and regular polygons", () => {
  test("roots of unity in radicals", () => {
    expect(String(rootOfUnityByRadicals(3))).toBe("(−1 + i·√3)/2");
    expect(String(rootOfUnityByRadicals(4))).toBe("i");
    expect(String(rootOfUnityByRadicals(5))).toBe("(−1 + √5 + √(−10 − 2·√5))/4");
    expect(String(rootOfUnityByRadicals(12))).toBe("(i + √3)/2");
    for (const n of [1, 2, 6, 7, 8, 9, 11, 13, 15, 16, 20]) expect(rootOfUnityByRadicals(n).check(30).agrees).toBe(true);
    expect(rootOfUnityByRadicals(7).indices).toEqual([2, 3]);
    expect(rootOfUnityByRadicals(11).indices).toEqual([2, 5]);
  });
  test("Gauss: the 17-gon in square roots", () => {
    const r = regularPolygon(17);
    expect(r.constructible).toBe(true);
    expect(String(r.cos)).toBe("(−1 + √17 + √(34 − 2·√17) + 2·√(17 + 3·√17 − √(170 + 38·√17)))/16");
    expect(r.cos.usesOnlySquareRoots).toBe(true);
    expect(r.cos.check(50).agrees).toBe(true);
    expect(r.cos.approx()).toBeCloseTo(Math.cos((2 * Math.PI) / 17), 14);
  });
  test("which polygons are constructible", () => {
    const yes = [];
    for (let n = 3; n <= 40; n += 1) if (regularPolygon(n).constructible) yes.push(n);
    expect(yes).toEqual([3, 4, 5, 6, 8, 10, 12, 15, 16, 17, 20, 24, 30, 32, 34, 40]);
    for (const n of yes) {
      const c = regularPolygon(n).cos;
      expect(c.usesOnlySquareRoots).toBe(true);
      expect(c.check(30).agrees).toBe(true);
      expect(c.approx()).toBeCloseTo(Math.cos((2 * Math.PI) / n), 13);
    }
    expect(String(regularPolygon(5).cos)).toBe("(−1 + √5)/4");
    expect(String(regularPolygon(16).cos)).toBe("√(2 + √2)/2");
    const seven = regularPolygon(7);
    expect(seven.constructible).toBe(false);
    expect(seven.cos.indices).toEqual([2, 3]);
    expect(regularPolygon(9).cos).toBe(null);
    expect(regularPolygon(257).cos.check(30).agrees).toBe(true);
  });
});

describe("ruler and compass", () => {
  test("constructible numbers", () => {
    const c = constructible("x^4 - 4x^2 + 2");
    expect(c.constructible).toBe(true);
    expect(c.factors[0].roots.map((r) => String(r.expression))).toEqual(["−√(2 + √2)", "−√(2 − √2)", "√(2 − √2)", "√(2 + √2)"]);
    expect(constructible("x^3 - 2").constructible).toBe(false);
    // degree 4, but the splitting field has degree 24
    const s4 = constructible("x^4 + x + 1");
    expect(s4.constructible).toBe(false);
    expect(s4.reason).toMatch(/24/);
    expect(constructible("x^8 - 2").constructible).toBe(null);
    expect(constructible("(x^2 - 2)*(x^3 - 2)").constructible).toBe("some");
  });
  test("the classical problems", () => {
    expect(classicalProblems.doublingTheCube().possible).toBe(false);
    const t = classicalProblems.trisecting();
    expect(t.possible).toBe(false);
    expect(t.polynomial).toBe("8x³ − 6x − 1");
    const ok = classicalProblems.trisecting("-11/16");
    expect(ok.possible).toBe(true);
    expect(ok.roots.map((r) => String(r.expression))).toEqual(["1/4", "(−1 − 3·√5)/8", "(−1 + 3·√5)/8"]);
    expect(classicalProblems.trisecting(0).possible).toBe(true);
    expect(classicalProblems.squaringTheCircle().possible).toBe(false);
    expect(() => classicalProblems.trisecting(2)).toThrow();
  });
});
