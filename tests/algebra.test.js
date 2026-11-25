import { describe, expect, test } from "bun:test";
import { Rational } from "@ratmath/core";
import {
  ring, gcd, lcm, content, primitivePart, squareFree, squareFreePart,
  det, sylvester, resultant, discriminant, jacobian, hessian,
  factor, isIrreducible, integerPrimitive, primeField, NumberField,
  elementary, powerSum, completeHomogeneous, isSymmetric, toElementary,
} from "../index.js";

const R = ring("x, y, z");
const P = (t) => R.poly(t);
const product = (fs) => fs.reduce((p, { factor: g, multiplicity }) => p.mul(g.pow(multiplicity)), R.one);

describe("gcd and lcm", () => {
  test("common factors in several variables", () => {
    const g = P("x^2 + y z - 1");
    const a = g.mul(P("x - y"));
    const b = g.mul(P("x + 2 z^2"));
    expect(gcd(a, b).equals(g)).toBe(true);
    expect(gcd([a, b, g.mul(3)]).equals(g)).toBe(true);
    expect(String(gcd(P("x^2 - y^2"), P("x^3 - y^3")))).toBe("x − y");
    expect(gcd(P("x + 1"), P("y + 1")).equals(R.one)).toBe(true);
  });
  test("zero, constants and numbers", () => {
    expect(gcd(R.zero, R.zero).isZero()).toBe(true);
    expect(gcd(R.zero, P("2x")).equals(P("x"))).toBe(true);
    expect(gcd(P("4x"), 6).equals(R.one)).toBe(true);
    expect(() => gcd(1, 2)).toThrow("polynomial");
    expect(() => lcm(1, 2)).toThrow("polynomial");
  });
  test("lcm times gcd is the product up to a constant", () => {
    const a = P("(x + y)(x - z)");
    const b = P("(x + y)(y + z)^2");
    expect(lcm(a, b).equals(P("(x + y)(x - z)(y + z)^2"))).toBe(true);
    expect(lcm(a, R.zero).isZero()).toBe(true);
    expect(lcm([a]).equals(a.monic())).toBe(true);
  });
  test("over a finite field and a number field", () => {
    const S = ring("x, y", { field: primeField(7) });
    const g = S.poly("x y + 3");
    expect(gcd(g.mul(S.poly("x + y")), g.mul(S.poly("x - y + 1"))).equals(g.monic())).toBe(true);
    const K = NumberField.quadratic(2);
    const T = ring("x", { field: K });
    const r = T.constant(K.gen);
    const x = T.gens[0];
    expect(gcd(T.poly("x^2 - 2"), x.sub(r).mul(T.poly("x + 5"))).equals(x.sub(r))).toBe(true);
  });
  test("content and primitive part", () => {
    const f = P("(y^2 - 1) x^2 + (y + 1) x + (y + 1) z");
    expect(String(content(f, "x"))).toBe("y + 1");
    expect(primitivePart(f, "x").mul(content(f, "x")).equals(f)).toBe(true);
    expect(content(R.zero, "x").isZero()).toBe(true);
    expect(primitivePart(R.zero, "x").isZero()).toBe(true);
  });
});

describe("square-free decomposition", () => {
  test("multiplicities are separated", () => {
    const f = P("3 (x + y) (x - z)^2 (y z + 1)^3");
    const { unit, factors } = squareFree(f);
    expect(unit).toEqual(new Rational(3));
    expect(factors.map((t) => t.multiplicity)).toEqual([1, 2, 3]);
    expect(product(factors).mul(3).equals(f)).toBe(true);
    expect(squareFreePart(f).equals(P("(x + y)(x - z)(y z + 1)").monic())).toBe(true);
  });
  test("factors of equal multiplicity stay together", () => {
    const { factors } = squareFree(P("(x^2 - 1)^2 y^2"));
    expect(factors).toHaveLength(1);
    expect(factors[0].factor.equals(P("x^2 y - y"))).toBe(true);
  });
  test("limits", () => {
    expect(() => squareFree(R.zero)).toThrow("zero polynomial");
    expect(() => squareFree(ring("x", { field: primeField(3) }).poly("x^3"))).toThrow("characteristic 0");
    expect(squareFree(P("5")).factors).toEqual([]);
  });
});

describe("determinants and resultants", () => {
  test("Vandermonde", () => {
    const [x, y, z] = R.gens;
    const d = det([[1, x, x.pow(2)], [1, y, y.pow(2)], [1, z, z.pow(2)]]);
    expect(d.equals(P("(y - x)(z - x)(z - y)"))).toBe(true);
  });
  test("pivoting, singular and trivial matrices", () => {
    expect(String(det([[0, P("x")], [P("y"), 0]]))).toBe("−xy");
    expect(det([[P("x"), P("y")], [P("2x"), P("2y")]]).isZero()).toBe(true);
    expect(det([[0, 0, 1], [0, 0, P("x")], [1, 2, 3]], R).isZero()).toBe(true);
    expect(det([], R).equals(R.one)).toBe(true);
    expect(det([[P("x + 1")]]).equals(P("x + 1"))).toBe(true);
    expect(() => det([[1, 2], [3, 4]])).toThrow("ring");
    expect(() => det([[P("x"), 1]])).toThrow("square");
    const [x] = R.gens;
    const m = [[x, 1, 0, 2], [1, x, 1, 0], [0, 1, x, 1], [2, 0, 1, x]];
    // against the permutation expansion
    const perms = (a) => (a.length <= 1 ? [[a, 1]] : a.flatMap((v, i) => perms([...a.slice(0, i), ...a.slice(i + 1)]).map(([p, s]) => [[v, ...p], i % 2 ? -s : s])));
    let leibniz = R.zero;
    for (const [p, s] of perms([0, 1, 2, 3])) leibniz = leibniz.add(p.reduce((t, j, i) => t.mul(m[i][j]), R.from(s)));
    expect(det(m).equals(leibniz)).toBe(true);
    expect(det(m).degree()).toBe(4);
  });
  test("the resultant eliminates a variable", () => {
    // circle and line: the x-coordinates of the intersections
    const r = resultant(P("x^2 + y^2 - 1"), P("x - y"), "y");
    expect(String(r)).toBe("2x² − 1");
    expect(resultant(P("x^2 - 1"), P("x - 1"), "x").isZero()).toBe(true);
    expect(resultant(P("x - 2"), P("x - 5"), "x").toConstant()).toEqual(new Rational(3));
    expect(resultant(P("y"), P("x^2 + 1"), "x").equals(P("y^2"))).toBe(true);
    expect(resultant(P("x^2 + 1"), P("y"), "x").equals(P("y^2"))).toBe(true);
    expect(resultant(R.zero, P("x"), "x").isZero()).toBe(true);
    expect(sylvester(P("x^2 + 1"), P("x - 3"), "x").map((row) => row.map(String))).toEqual([["1", "0", "1"], ["1", "−3", "0"], ["0", "1", "−3"]]);
    expect(() => sylvester(R.zero, P("x"), "x")).toThrow("zero polynomial");
  });
  test("discriminants", () => {
    const S = ring("a, b, c, x");
    expect(String(discriminant(S.poly("a x^2 + b x + c"), "x"))).toBe("b² − 4ac");
    const T = ring("p, q, x");
    expect(String(discriminant(T.poly("x^3 + p x + q"), "x"))).toBe("−4p³ − 27q²");
    expect(() => discriminant(P("y"), "x")).toThrow("degree at least 1");
  });
  test("jacobian and hessian", () => {
    const S = ring("x, y");
    const J = jacobian([S.poly("x^2 y"), S.poly("x + y^3")]);
    expect(J.map((row) => row.map(String))).toEqual([["2xy", "x²"], ["1", "3y²"]]);
    expect(String(det(J))).toBe("6xy³ − x²");
    expect(hessian(S.poly("x^3 + x y^2")).map((row) => row.map(String))).toEqual([["6x", "2y"], ["2y", "2x"]]);
  });
});

describe("factoring over Q", () => {
  const check = (text, expected) => {
    const f = P(text);
    const { unit, factors } = factor(f);
    expect(product(factors).mul(unit).equals(f)).toBe(true);
    expect(factors.map((t) => `${t.factor}${t.multiplicity > 1 ? `^${t.multiplicity}` : ""}`)).toEqual(expected);
  };
  test("classic identities", () => {
    check("x^2 - y^2", ["x + y", "x − y"]);
    check("x^3 + y^3 + z^3 - 3 x y z", ["x + y + z", "x² − xy + y² − xz − yz + z²"]);
    check("x^4 + 4 y^4", ["x² + 2xy + 2y²", "x² − 2xy + 2y²"]);
  });
  test("units, contents and multiplicities", () => {
    const f = P("-3/2 (x + y)^2 (x y - 1) z");
    const { unit, factors } = factor(f);
    expect(unit).toEqual(new Rational(-3, 2));
    expect(product(factors).mul(unit).equals(f)).toBe(true);
    expect(factors.map((t) => t.multiplicity).sort()).toEqual([1, 1, 2]);
    expect(factor(P("7")).factors).toEqual([]);
    expect(factor(P("7")).unit).toEqual(new Rational(7));
  });
  test("one variable goes through the same door", () => {
    check("x^4 - 1", ["x + 1", "x − 1", "x² + 1"]);
    check("2 y^2 + 3 y + 1", ["2y + 1", "y + 1"]);
  });
  test("irreducibility", () => {
    expect(isIrreducible(P("x^2 + y^2 - 1"))).toBe(true);
    expect(isIrreducible(P("x^2 + y^2"))).toBe(true);
    expect(isIrreducible(P("x y - z^2"))).toBe(true);
    expect(isIrreducible(P("x^2 - y^2"))).toBe(false);
    expect(isIrreducible(P("(x + y)^2"))).toBe(false);
    expect(isIrreducible(P("3"))).toBe(false);
  });
  test("integerPrimitive", () => {
    const { content: c, primitive } = integerPrimitive(P("-2/3 x + 4/9 y"));
    expect(c).toEqual(new Rational(-2, 9));
    expect(String(primitive)).toBe("3x − 2y");
  });
  test("limits", () => {
    expect(() => factor(R.zero)).toThrow("zero polynomial");
    expect(() => factor(ring("x", { field: primeField(5) }).poly("x^2 - 1"))).toThrow("rationals");
  });
});

describe("symmetric polynomials", () => {
  test("the standard families", () => {
    expect(String(elementary(R, 2))).toBe("xy + xz + yz");
    expect(elementary(R, 0).equals(R.one)).toBe(true);
    expect(elementary(R, 4).isZero()).toBe(true);
    expect(String(powerSum(R, 3))).toBe("x³ + y³ + z³");
    expect(completeHomogeneous(R, 2).equals(P("x^2 + y^2 + z^2 + x y + x z + y z"))).toBe(true);
    expect(completeHomogeneous(R, 0).equals(R.one)).toBe(true);
  });
  test("isSymmetric", () => {
    expect(isSymmetric(P("x y z + x + y + z"))).toBe(true);
    expect(isSymmetric(P("x^2 y + y^2 z + z^2 x"))).toBe(false); // cyclic only
    expect(isSymmetric(P("x + y"))).toBe(false);
    expect(isSymmetric(ring("x").poly("x^3 + 1"))).toBe(true);
  });
  test("Newton's identities fall out of toElementary", () => {
    expect(String(toElementary(powerSum(R, 2)).polynomial)).toBe("e1² − 2e2");
    expect(String(toElementary(powerSum(R, 3)).polynomial)).toBe("e1³ − 3e1e2 + 3e3".replace("e1e2", "e1·e2"));
  });
  test("substituting back recovers the polynomial", () => {
    const f = P("(x - y)^2 (x - z)^2 (y - z)^2");
    const { polynomial, ring: E } = toElementary(f, { prefix: "s" });
    expect(E.vars).toEqual(["s1", "s2", "s3"]);
    const back = polynomial.compose({ s1: elementary(R, 1), s2: elementary(R, 2), s3: elementary(R, 3) });
    expect(back.equals(f)).toBe(true);
    expect(() => toElementary(P("x - y"))).toThrow("not symmetric");
    expect(toElementary(R.zero).polynomial.isZero()).toBe(true);
  });
});
