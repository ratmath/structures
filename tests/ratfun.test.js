import { describe, expect, test } from "bun:test";
import { Rational } from "@ratmath/core";
import { ring, fractions, RationalFunction, primeField, det as polyDet } from "../index.js";

const R = ring("x, y");
const F = fractions(R);
const f = (t) => F.parse(t);

// the same Gauss–Jordan elimination @ratmath/linalg runs over any field object
function inverse(K, m) {
  const n = m.length;
  const a = m.map((row, i) => [...row, ...Array.from({ length: n }, (_, j) => (i === j ? K.one : K.zero))]);
  for (let c = 0; c < n; c += 1) {
    const p = a.findIndex((row, i) => i >= c && K.sign(row[c]) !== 0);
    if (p < 0) return null;
    [a[c], a[p]] = [a[p], a[c]];
    const lead = a[c][c];
    a[c] = a[c].map((v) => K.div(v, lead));
    for (let i = 0; i < n; i += 1) {
      if (i === c) continue;
      const k = a[i][c];
      a[i] = a[i].map((v, j) => K.sub(v, K.mul(k, a[c][j])));
    }
  }
  return a.map((row) => row.slice(n));
}

describe("rational functions", () => {
  test("lowest terms and a canonical form", () => {
    expect(String(f("(x^2 - y^2)/(x - y)"))).toBe("x + y");
    expect(String(f("(x^2 - 1)/(x^2 + 2x + 1)"))).toBe("(x − 1)/(x + 1)");
    expect(String(f("1/(2x) + 1/(3y)"))).toBe("(2x + 3y)/(6xy)");
    expect(String(f("x/(-y)"))).toBe("−x/y");
    expect(f("(1/2 x)/(1/3 y)").equals(f("3x/(2y)"))).toBe(true);
    expect(F.name).toBe("Q(x, y)");
    expect(String(F)).toBe("Q(x, y)");
    expect(fractions(R)).toBe(F);
  });
  test("field operations", () => {
    const a = f("1/(x + y)");
    const b = f("1/(x - y)");
    expect(a.add(b).equals(f("2x/(x^2 - y^2)"))).toBe(true);
    expect(a.sub(b).equals(f("-2y/(x^2 - y^2)"))).toBe(true);
    expect(a.mul(b).equals(f("1/(x^2 - y^2)"))).toBe(true);
    expect(a.div(b).equals(f("(x - y)/(x + y)"))).toBe(true);
    expect(a.inv().isPolynomial()).toBe(true);
    expect(a.inv().toPolynomial().equals(R.poly("x + y"))).toBe(true);
    expect(() => a.toPolynomial()).toThrow();
    expect(a.pow(-2).equals(f("(x + y)^2"))).toBe(true);
    expect(a.pow(0).equals(F.one)).toBe(true);
    expect(a.neg().add(a).isZero()).toBe(true);
    expect(() => F.zero.inv()).toThrow("division by zero");
    expect(() => a.div(0)).toThrow("division by zero");
    expect(() => F.of(1, 0)).toThrow("division by zero");
    expect(a.add(1).equals(f("(x + y + 1)/(x + y)"))).toBe(true);
    expect(a.mul("x + y").equals(F.one)).toBe(true);
  });
  test("negative exponents in text", () => {
    expect(f("x^-2 y").equals(f("y/x^2"))).toBe(true);
    expect(f("(x/y)^-1").equals(f("y/x"))).toBe(true);
    expect(() => f("x^y")).toThrow("exponent");
    expect(() => f("1/0")).toThrow("division by zero");
  });
  test("derivatives, substitution and values", () => {
    expect(f("x/(x + y)").diff("x").equals(f("y/(x + y)^2"))).toBe(true);
    expect(f("1/x").diff("y").isZero()).toBe(true);
    expect(f("(x + y)/(x - y)").subs({ y: 1 }).equals(f("(x + 1)/(x - 1)"))).toBe(true);
    expect(f("(x + y)/(x - y)").at({ x: 3, y: 1 })).toEqual(new Rational(2));
    expect(() => f("1/(x - y)").at({ x: 1, y: 1 })).toThrow();
    expect(() => f("1/(x - y)").subs({ x: "y" })).toThrow();
  });
  test("printing", () => {
    expect(f("(x + 1)/y").toText()).toBe("(x + 1)/y");
    expect(F.parse(f("(x + 1)/(y^2 - 2)").toText()).equals(f("(x + 1)/(y^2 - 2)"))).toBe(true);
    expect(f("(x + 1)/(2 y)").toTeX()).toBe("\\frac{x + 1}{2y}");
    expect(f("x + 1").toTeX()).toBe("x + 1");
    expect(F.format(f("1/x"))).toBe("1/x");
    expect(String(f("3/(x y)"))).toBe("3/(xy)");
    expect(String(f("x/y^2"))).toBe("x/(y²)");
    expect(String(F.zero)).toBe("0");
    // a constant denominator is still a polynomial
    expect(String(f("x/2 + 1/3"))).toBe("(1/2)x + 1/3");
    expect(f("x/2 + 1/3").toText()).toBe("1/2*x + 1/3");
    expect(f("x/2").toTeX()).toBe("\\frac{1}{2}x");
    expect(f("(x + y)/2").toPolynomial().equals(R.poly("1/2 x + 1/2 y"))).toBe(true);
    expect(f("(x + y)/2").at({ x: 1, y: 2 })).toEqual(new Rational(3, 2));
  });
  test("from", () => {
    expect(F.from(R.poly("x")).equals(F.gens[0])).toBe(true);
    expect(F.from("1/2").equals(F.of(1, 2))).toBe(true);
    expect(F.from(3).equals(f("3"))).toBe(true);
    expect(F.from(F.one)).toBe(F.one);
    const G = fractions(ring("x, y, t"));
    expect(String(G.from(f("x/y")))).toBe("x/y");
    expect(f("x/y") instanceof RationalFunction).toBe(true);
    expect(f("x/y").ring).toBe(R);
  });
  test("other coefficient fields keep the denominator monic", () => {
    const S = ring("t", { field: primeField(5) });
    const G = fractions(S);
    expect(String(G.parse("(t + 1)/(2t + 1)"))).toBe("(3t + 3)/(t + 3)");
    expect(G.parse("1/t + 1/t + 1/t + 1/t + 1/t").isZero()).toBe(true);
    expect(G.name).toBe("F_5(t)");
  });
});

describe("matrices with parameters", () => {
  test("inverse of a symbolic 2×2 matrix", () => {
    const S = ring("a, b, c, d");
    const K = fractions(S);
    const [a, b, c, d] = K.gens;
    const inv = inverse(K, [[a, b], [c, d]]);
    const D = K.parse("a d - b c");
    expect(inv[0][0].equals(d.div(D))).toBe(true);
    expect(inv[0][1].equals(b.neg().div(D))).toBe(true);
    expect(inv[1][0].equals(c.neg().div(D))).toBe(true);
    expect(inv[1][1].equals(a.div(D))).toBe(true);
  });
  test("a singular family is detected", () => {
    const K = fractions(ring("t"));
    expect(inverse(K, [[K.parse("t"), K.parse("t^2")], [K.one, K.parse("t")]])).toBe(null);
  });
  test("characteristic polynomial through the fraction-free determinant", () => {
    const S = ring("t");
    const A = [[2, 1, 0], [1, 2, 1], [0, 1, 2]];
    const M = A.map((row, i) => row.map((v, j) => (i === j ? S.poly("t").sub(v) : S.from(-v))));
    expect(String(polyDet(M))).toBe("t³ − 6t² + 10t − 4");
  });
});
