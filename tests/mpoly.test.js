import { describe, expect, test } from "bun:test";
import { Rational } from "@ratmath/core";
import { ring, QQ, primeField, galoisField, NumberField, PolyRing } from "../index.js";

const R = ring("x, y, z");
const P = (t) => R.poly(t);

describe("rings", () => {
  test("variables from text or a list", () => {
    expect(ring("x y").vars).toEqual(["x", "y"]);
    expect(ring(["a", "b1"]).vars).toEqual(["a", "b1"]);
    expect(String(R)).toBe("Q[x, y, z]");
    expect(String(ring("x", { field: primeField(5) }))).toBe("F_5[x]");
    expect(String(ring("x", { field: NumberField.quadratic(2) }))).toBe("Q(√2)[x]");
  });
  test("bad variable lists are refused", () => {
    expect(() => ring("")).toThrow();
    expect(() => ring("x, x")).toThrow("distinct");
    expect(() => ring("2x")).toThrow("variable name");
    expect(() => ring("x", { order: "nope" })).toThrow("monomial order");
  });
  test("generators and constants", () => {
    const [x, y] = R.gens;
    expect(String(x.mul(y).add(1))).toBe("xy + 1");
    expect(R.gen("z")).toBe(R.gens[2]);
    expect(() => R.gen("w")).toThrow("not a variable");
    expect(R.from(0).isZero()).toBe(true);
    expect(R.from("3/4").toConstant()).toEqual(new Rational(3, 4));
    expect(R.from(new Rational(1, 2)).mul(4).toConstant()).toEqual(new Rational(2));
  });
});

describe("reading and printing", () => {
  test("juxtaposition, powers and signs", () => {
    expect(String(P("x^2 y - 3/2 x z + 1"))).toBe("x²y − (3/2)xz + 1");
    expect(String(P("3xy^2"))).toBe("3xy²");
    expect(String(P("-x^2"))).toBe("−x²");
    expect(String(P("(x+1)(x-1)"))).toBe("x² − 1");
    expect(String(P("2(x + y)^2"))).toBe("2x² + 4xy + 2y²");
    expect(String(P("x**3 − 2·x × y"))).toBe("x³ − 2xy");
    expect(String(P("- - x + + y"))).toBe("x + y");
    expect(String(P("x * -y"))).toBe("−xy");
    expect(String(P("x^(1+1)"))).toBe("x²");
  });
  test("exact numbers", () => {
    expect(P("1.25 x").coefficient({ x: 1 })).toEqual(new Rational(5, 4));
    expect(P("0.#3").toConstant()).toEqual(new Rational(1, 3));
    expect(P("x/3 + 2/4").equals(P("1/3 x + 1/2"))).toBe(true);
    expect(P("6 x y / (2 y)").equals(P("3x"))).toBe(true);
  });
  test("errors say where", () => {
    expect(() => P("x +")).toThrow("ends too soon");
    expect(() => P("x ^ y")).toThrow("exponent");
    expect(() => P("x ^ -1")).toThrow("exponent");
    expect(() => P("w + 1")).toThrow('"w" is not a variable');
    expect(() => P("x / y")).toThrow("does not divide");
    expect(() => P("x / 0")).toThrow("division by zero");
    expect(() => P("(x")).toThrow('expected ")"');
    expect(() => P("x $ y")).toThrow("unexpected");
    expect(() => P("x^2^3")).toThrow("parentheses");
    expect(() => P("")).toThrow("nothing to read");
  });
  test("text round trip and TeX", () => {
    const f = P("x^2 y - 3/2 x z + 1");
    expect(f.toText()).toBe("x^2*y - 3/2*x*z + 1");
    expect(P(f.toText()).equals(f)).toBe(true);
    expect(f.toTeX()).toBe("x^{2}y - \\frac{3}{2}xz + 1");
    expect(String(R.zero)).toBe("0");
    expect(ring("x1 x2").poly("x1^2 x2 + 2 x2").toString()).toBe("x1²·x2 + 2x2");
    expect(ring("x1 x2").poly("x1^2 x2").toTeX()).toBe("x_{1}^{2}\\,x_{2}");
  });
});

describe("arithmetic", () => {
  test("ring laws on examples", () => {
    const f = P("x^2 - y z + 3"), g = P("x + y - 1/2"), h = P("z^3 - x");
    expect(f.mul(g.add(h)).equals(f.mul(g).add(f.mul(h)))).toBe(true);
    expect(f.mul(g).equals(g.mul(f))).toBe(true);
    expect(f.sub(f).isZero()).toBe(true);
    expect(f.add(f.neg()).length).toBe(0);
    expect(g.pow(3).equals(g.mul(g).mul(g))).toBe(true);
    expect(g.pow(0).equals(R.one)).toBe(true);
    expect(() => g.pow(-1)).toThrow();
  });
  test("mixed operands", () => {
    expect(String(P("x").add(1).mul("y - 1"))).toBe("xy − x + y − 1");
    expect(String(P("x").scale("2/3"))).toBe("(2/3)x");
  });
  test("degrees and leading terms depend on the order", () => {
    const f = "x y^2 + x^2 + y^3 z";
    expect(String(R.poly(f).leadingTerm())).toBe("y³z");
    expect(String(R.withOrder("lex").poly(f).leadingTerm())).toBe("x²");
    expect(String(R.withOrder("grlex").poly("x y^2 z + x^2 z^2").leadingTerm())).toBe("x²z²");
    expect(String(R.withOrder("grevlex").poly("x y^2 z + x^2 z^2").leadingTerm())).toBe("xy²z");
    const g = P(f);
    expect(g.degree()).toBe(4);
    expect(g.degree("x")).toBe(2);
    expect(g.degrees()).toEqual({ x: 2, y: 3, z: 1 });
    expect(R.zero.degree()).toBe(-Infinity);
    expect(g.variables()).toEqual(["x", "y", "z"]);
    expect(P("x^2 + 1").variables()).toEqual(["x"]);
    expect(g.leadingExponents()).toEqual([0, 3, 1]);
    expect(() => R.zero.leadingTerm()).toThrow();
  });
  test("a custom order", () => {
    const S = ring("x y", { order: (a, b) => a[1] - b[1] || a[0] - b[0] });
    expect(String(S.poly("x^5 + y").leadingTerm())).toBe("y");
  });
  test("terms and coefficients", () => {
    const f = P("2 x^2 y - 3 y + 5");
    expect(f.terms().map((t) => [String(t.coefficient), t.exponents])).toEqual([["2", [2, 1, 0]], ["-3", [0, 1, 0]], ["5", [0, 0, 0]]]);
    expect(f.coefficient([2, 1, 0])).toEqual(new Rational(2));
    expect(f.coefficient({ y: 1 })).toEqual(new Rational(-3));
    expect(f.coefficient({ z: 4 })).toEqual(new Rational(0));
    expect(f.constantTerm()).toEqual(new Rational(5));
    expect(f.coefficients("y").map(String)).toEqual(["5", "2x² − 3"]);
    expect(f.isHomogeneous()).toBe(false);
    expect(P("x^2 + y z").isHomogeneous()).toBe(true);
    expect(String(P("x^2 + y + 1").homogenize())).toBe("x² + yh + h²");
    expect(R.fromTerms([[2, { x: 1 }], ["1/2", [0, 0, 3]]]).equals(P("2x + 1/2 z^3"))).toBe(true);
  });
  test("univariate round trip", () => {
    expect(P("3 y^2 - 1").toUnivariate().map(String)).toEqual(["-1", "0", "3"]);
    expect(R.fromUnivariate([1, 0, "1/2"], "z").equals(P("1/2 z^2 + 1"))).toBe(true);
    expect(() => P("x y").toUnivariate()).toThrow("more than one");
    expect(R.zero.toUnivariate()).toEqual([]);
  });
});

describe("division", () => {
  test("one divisor", () => {
    const { quotient, remainder } = P("x^3 - y^3").divmod(P("x - y"));
    expect(String(quotient)).toBe("x² + xy + y²");
    expect(remainder.isZero()).toBe(true);
    expect(P("x^3 - y^3").isDivisibleBy("x - y")).toBe(true);
    expect(P("x^3 - y^3").isDivisibleBy("x + y")).toBe(false);
    expect(() => P("x^2 + 1").div("x - y")).toThrow("does not divide");
    expect(() => P("x").divmod(R.zero)).toThrow("zero polynomial");
    expect(String(P("4x").div(2))).toBe("2x");
  });
  test("several divisors, the textbook example", () => {
    // Cox, Little, O'Shea, chapter 2 §3, with lex order x > y
    const L = ring("x, y", { order: "lex" });
    const f = L.poly("x^2 y + x y^2 + y^2");
    const { quotients, remainder } = f.divmod([L.poly("x y - 1"), L.poly("y^2 - 1")]);
    expect(quotients.map(String)).toEqual(["x + y", "1"]);
    expect(String(remainder)).toBe("x + y + 1");
    // the other order of divisors gives another remainder
    const other = f.divmod([L.poly("y^2 - 1"), L.poly("x y - 1")]);
    expect(String(other.remainder)).toBe("2x + 1");
    for (const r of [{ quotients, remainder, gs: ["x y - 1", "y^2 - 1"] }, { ...other, gs: ["y^2 - 1", "x y - 1"] }]) {
      const back = r.quotients.reduce((s, q, i) => s.add(q.mul(L.poly(r.gs[i]))), r.remainder);
      expect(back.equals(f)).toBe(true);
    }
  });
});

describe("calculus and substitution", () => {
  test("partial derivatives", () => {
    const f = P("x^3 y^2 - 2 x z + 7");
    expect(String(f.diff("x"))).toBe("3x²y² − 2z");
    expect(String(f.diff("x", 2))).toBe("6xy²");
    expect(f.diff("x").diff("y").equals(f.diff("y").diff("x"))).toBe(true);
    expect(f.gradient().map(String)).toEqual(["3x²y² − 2z", "2x³y", "−2x"]);
  });
  test("substitution and evaluation", () => {
    const f = P("x^2 y - 3/2 x z + 1");
    expect(String(f.subs({ x: 2 }))).toBe("4y − 3z + 1");
    expect(String(f.subs({ x: "y + 1", z: 0 }))).toBe("y³ + 2y² + y + 1");
    expect(f.at({ x: 2, y: "1/3", z: 1 })).toEqual(new Rational(-2, 3));
    expect(f.at([1, 1, 1])).toEqual(new Rational(1, 2));
    expect(() => f.at({ x: 1 })).toThrow("no value was given for y");
    expect(P("x + 1").at({ x: 4 })).toEqual(new Rational(5));
    // simultaneous, not one after another
    expect(String(P("x - y").subs({ x: "y", y: "x" }))).toBe("−x + y");
  });
  test("moving between rings", () => {
    const S = ring("y, x, w", { order: "lex" });
    expect(String(P("x^2 y + 1").into(S))).toBe("yx² + 1");
    expect(() => P("z").into(S)).toThrow("no variable z");
    const F = ring("x y", { field: primeField(3) });
    expect(String(P("4 x y + 3 x - 1/2").into(F))).toBe("xy + 1");
  });
});

describe("other fields", () => {
  test("integers mod p", () => {
    const F = ring("x y", { field: primeField(5) });
    expect(String(F.poly("(x + y)^5"))).toBe("x⁵ + y⁵");
    expect(String(F.poly("3x - 7/2 y"))).toBe("3x + 4y");
    expect(F.poly("x^5").diff("x").isZero()).toBe(true);
    expect(F.poly("2 x").at({ x: 4 })).toBe(3n);
  });
  test("a field with p^k elements", () => {
    const G = galoisField(2, 2);
    const F = ring("x", { field: G });
    const a = F.constant(G.gen);
    // x² + x + 1 splits over GF(4)
    const f = F.poly("x").add(a).mul(F.poly("x").add(a).add(1));
    expect(f.equals(F.poly("x^2 + x + 1"))).toBe(true);
  });
  test("a number field", () => {
    const K = NumberField.quadratic(2);
    const S = ring("x y", { field: K });
    const r = S.constant(K.gen);
    expect(String(S.poly("x").sub(r).mul(S.poly("x").add(r)))).toBe("x² − 2");
    expect(String(S.poly("x + y").mul(r))).toBe("√2·x + √2·y");
    expect(String(S.poly("x").add(r).pow(2))).toBe("x² + 2√2·x + 2");
  });
  test("a field must be exact", () => {
    const fuzzy = { zero: 0, one: 1, add: (a, b) => a + b, sub: (a, b) => a - b, mul: (a, b) => a * b, div: (a, b) => a / b, neg: (a) => -a, sign: () => null };
    expect(() => ring("x", { field: fuzzy }).poly("x - x")).toThrow("cannot decide");
    expect(() => new PolyRing("x", { field: {} })).toThrow();
    expect(ring("x", { field: QQ }).K).toBe(QQ);
  });
});
