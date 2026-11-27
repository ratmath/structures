import { describe, expect, test } from "bun:test";
import { extension, conwayPolynomial, conwayField, formatModular, identify } from "../index.js";

describe("towers of fields", () => {
  const K = extension({ a: "a^2 - 2", b: "b^2 - 3" });
  test("degrees and the tower law", () => {
    expect(K.degree).toBe(4);
    expect(K.degrees).toEqual([2, 2]);
    expect(K.towerLaw()).toBe("[Q(a, b) : Q] = [Q(a, b) : Q(a)]·[Q(a) : Q] = 2·2 = 4");
    expect(K.basis.map(String)).toEqual(["1", "b", "a", "ab"]);
    const M = extension({ a: "a^2 - 2", b: "b^3 - a" });
    expect(M.degree).toBe(6);
    expect(M.degreeOver(1)).toBe(3);
    expect(M.b.minimalPolynomialText()).toBe("x⁶ − 2");
    expect(M.of("b^3").equals(M.a)).toBe(true);
  });
  test("arithmetic in the tower basis", () => {
    const x = K.of("a + b");
    expect(String(x.mul(x))).toBe("2ab + 5");
    expect(String(x.inverse())).toBe("−a + b");
    expect(String(K.of("1/(a + b)"))).toBe("−a + b");
    expect(x.div(x).equals(1)).toBe(true);
    expect(x.pow(4).sub(x.pow(2).mul(10)).add(1).isZero()).toBe(true);
    expect(x.minimalPolynomialText()).toBe("x⁴ − 10x² + 1");
    expect(x.degree).toBe(4);
    expect(K.of("a*b").degree).toBe(2);
    expect(String(x.norm())).toBe("1");
    expect(String(x.trace())).toBe("0");
    expect(x.approx()).toBeCloseTo(Math.SQRT2 + Math.sqrt(3), 12);
    expect(K.a.approx()).toBeCloseTo(Math.SQRT2, 12);
    expect(x.coordinates().map(String)).toEqual(["0", "1", "1", "0"]);
  });
  test("a reducible step is refused, with the reason", () => {
    expect(() => extension({ a: "a^2 - 2", b: "b^2 - 8" })).toThrow(/not irreducible over Q\(a\)/);
    expect(() => extension({ a: "a^2 - 4" })).toThrow(/not irreducible over Q/);
    expect(() => extension({ a: "a^2 - b", b: "b^2 - 2" })).toThrow(/comes later/);
  });
  test("primitive elements, roots and factoring over the field", () => {
    const p = K.primitiveElement();
    expect(p.text).toBe("x⁴ − 10x² + 1");
    expect(String(p.element)).toBe("a + b");
    expect(K.roots("x^2 - 6").map((r) => String(r.value)).sort()).toEqual(["ab", "−ab"]);
    expect(K.roots("x^2 - 5")).toEqual([]);
    expect(K.factor("x^4 - 10x^2 + 1").map((f) => f.degree)).toEqual([1, 1, 1, 1]);
    expect(K.factor("x^4 - 2").map((f) => f.text).sort()).toEqual(["x² + a", "x² − a"]);
    expect(K.subfield("a*b").degree).toBe(2);
    expect(K.subfield("a + b").isWhole).toBe(true);
  });
  test("automorphisms and the Galois correspondence", () => {
    expect(K.isGalois).toBe(true);
    const G = K.galoisGroup();
    expect(G.order).toBe(4);
    expect(identify(G).name).toBe("C2 × C2");
    const fixed = K.automorphisms().slice(1).map((s) => K.fixedField(s));
    expect(fixed.map((F) => F.degree)).toEqual([2, 2, 2]);
    expect(fixed.map((F) => F.text).sort()).toEqual(["x² − 2", "x² − 3", "x² − 6"]);
    expect(K.fixedField(K.automorphisms()).degree).toBe(1);

    const S = extension({ c: "c^3 - 2", w: "w^2 + w + 1" });
    expect(S.isGalois).toBe(true);
    const GS = S.galoisGroup();
    expect(GS.order).toBe(6);
    expect(GS.isAbelian).toBe(false);
    // the subgroup of order 3 fixes Q(w); the three of order 2 fix the three cubic fields
    const degrees = S.automorphisms().map((s) => S.fixedField(s).degree);
    expect(degrees.sort()).toEqual([2, 2, 3, 3, 3, 6]);
    const s = S.automorphisms()[1];
    expect(s.apply(S.c.mul(S.w)).equals(s.apply(S.c).mul(s.apply(S.w)))).toBe(true);

    const M = extension({ c: "c^3 - 2" });
    expect(M.isGalois).toBe(false);
    expect(M.automorphisms().length).toBe(1);
    expect(() => M.galoisGroup()).toThrow(/not a Galois extension/);
  });
  test("minimal polynomials over a field in between", () => {
    const M = extension({ a: "a^2 - 2", b: "b^3 - a" });
    const m = M.of("b").minimalPolynomialOver(1);
    expect(m.map(String)).toEqual(["−a", "0", "0", "1"]);
    expect(M.of("b").minimalPolynomialOver(0).length).toBe(7);
    expect(M.of("a + 1").minimalPolynomialOver(1).map(String)).toEqual(["−a − 1", "1"]);
  });
});

describe("Conway polynomials", () => {
  // values from Frank Lübeck's published tables
  const known = {
    "2,1": "x + 1", "2,2": "x² + x + 1", "2,3": "x³ + x + 1", "2,4": "x⁴ + x + 1", "2,5": "x⁵ + x² + 1", "2,6": "x⁶ + x⁴ + x³ + x + 1",
    "2,7": "x⁷ + x + 1", "2,8": "x⁸ + x⁴ + x³ + x² + 1", "3,1": "x + 1", "3,2": "x² + 2x + 2", "3,3": "x³ + 2x + 1", "3,4": "x⁴ + 2x³ + 2",
    "3,5": "x⁵ + 2x + 1", "3,6": "x⁶ + 2x⁴ + x² + 2x + 2", "5,1": "x + 3", "5,2": "x² + 4x + 2", "5,3": "x³ + 3x + 3", "5,4": "x⁴ + 4x² + 4x + 2",
    "7,1": "x + 4", "7,2": "x² + 6x + 3", "7,3": "x³ + 6x² + 4", "11,2": "x² + 7x + 2", "13,2": "x² + 12x + 2",
  };
  test("the published values are found by search", () => {
    for (const [k, v] of Object.entries(known)) {
      const [p, n] = k.split(",").map(Number);
      expect(formatModular(conwayPolynomial(p, n))).toBe(v);
    }
  });
  test("the standard generator is primitive and the fields fit together", () => {
    const F = conwayField(2, 6);
    expect(F.order(F.gen)).toBe(63n);
    for (const m of [1, 2, 3]) {
      const e = F.embed(m);
      const small = e.from;
      // the image of the small generator satisfies the small Conway polynomial
      let acc = F.zero;
      for (let i = small.conway.length - 1; i >= 0; i -= 1) acc = F.add(F.mul(acc, e.image), F.from(small.conway[i]));
      expect(F.sign(acc)).toBe(0);
    }
    const e = F.embed(3);
    const s = e.from;
    const x = [1n, 1n, 0n];
    const y = [0n, 1n, 1n];
    expect(F.equals(e.map(s.mul(x, y)), F.mul(e.map(x), e.map(y)))).toBe(true);
    expect(() => F.embed(4)).toThrow();
    expect(() => conwayPolynomial(4, 2)).toThrow();
    expect(() => conwayPolynomial(2, 40)).toThrow(/limit/);
  });
});
