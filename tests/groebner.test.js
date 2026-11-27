import { describe, expect, test } from "bun:test";
import { ring, groebner, sPolynomial, eliminationOrder, ideal, solve, chineseRemainder } from "../index.js";

const R = ring("x, y");
const strs = (xs) => xs.map(String);

const cyclic = (n) => {
  const S = ring(Array.from({ length: n }, (_, i) => `x${i}`));
  const g = S.gens;
  const out = [];
  for (let k = 1; k < n; k += 1) {
    let s = S.zero;
    for (let i = 0; i < n; i += 1) {
      let m = S.one;
      for (let j = 0; j < k; j += 1) m = m.mul(g[(i + j) % n]);
      s = s.add(m);
    }
    out.push(s);
  }
  out.push(g.reduce((m, v) => m.mul(v), S.one).sub(S.one));
  return out;
};

describe("Gröbner bases", () => {
  test("a small lex basis is triangular", () => {
    const S = ring("x, y, z", { order: "lex" });
    const G = groebner([S.poly("x^2 + y^2 + z^2 - 1"), S.poly("x - y"), S.poly("y - z^2")]);
    expect(strs(G.basis)).toEqual(["z⁴ + (1/2)z² − 1/2", "y − z²", "x − z²"]);
    expect(G.verify().ok).toBe(true);
    expect(G.dimension).toBe(0);
    expect(G.quotientDimension).toBe(4);
  });
  test("membership, remainders and certificates", () => {
    const f1 = R.poly("x^2 + y^2 - 1");
    const f2 = R.poly("x*y - 1/2");
    const G = groebner([f1, f2]);
    expect(strs(G.basis)).toEqual(["xy − 1/2", "x² + y² − 1", "y³ + (1/2)x − y"]);
    expect(strs(G.standardMonomials())).toEqual(["1", "y", "x", "y²"]);
    const f = f1.mul(R.poly("x + 3")).add(f2.mul(R.poly("y^2")));
    expect(G.contains(f)).toBe(true);
    expect(G.reduce(f).isZero()).toBe(true);
    expect(G.contains(R.poly("x"))).toBe(false);
    const ex = G.express(f);
    expect(ex.check).toBe(true);
    expect(ex.cofactors[0].mul(f1).add(ex.cofactors[1].mul(f2)).equals(f)).toBe(true);
    expect(() => G.express(R.poly("x"))).toThrow();
  });
  test("the trivial and the zero ideal", () => {
    expect(groebner([R.poly("x"), R.poly("x - 1")]).isTrivial).toBe(true);
    expect(groebner([R.zero]).isZero).toBe(true);
    expect(groebner([R.poly("x*y")]).dimension).toBe(1);
  });
  test("S-polynomials cancel the leading terms", () => {
    const s = sPolynomial(R.poly("x^2 y - 1"), R.poly("x y^2 - x"));
    expect(String(s)).toBe("x² − y");
  });
  test("an elimination order removes variables", () => {
    const S = ring("t, x, y");
    const G = groebner([S.poly("x - t^2"), S.poly("y - t^3")], { order: eliminationOrder(3, [0]) });
    const free = G.basis.filter((g) => g.degree("t") === 0);
    expect(strs(free)).toEqual(["x³ − y²"]);
  });
  test("benchmarks: cyclic 4, cyclic 5", () => {
    const G4 = groebner(cyclic(4));
    expect(G4.length).toBe(7);
    expect(G4.dimension).toBe(1);
    const G5 = groebner(cyclic(5));
    expect(G5.length).toBe(20);
    expect(G5.quotientDimension).toBe(70);
    expect(G5.verify().ok).toBe(true);
  });
  test("over a finite field", async () => {
    const { primeField } = await import("../index.js");
    const S = ring("x, y", { field: primeField(5) });
    const G = groebner([S.poly("x^2 + y^2 - 1"), S.poly("x y - 2")]);
    expect(G.verify().ok).toBe(true);
    expect(G.quotientDimension).toBe(4);
  });
});

describe("ideals", () => {
  const I = ideal(R, ["x"]);
  const J = ideal(R, ["y"]);
  test("sum, product, intersection, quotient, saturation", () => {
    expect(strs(I.intersect(J).basis)).toEqual(["xy"]);
    expect(strs(I.add(J).basis)).toEqual(["y", "x"]);
    expect(strs(I.mul(J).basis)).toEqual(["xy"]);
    expect(strs(ideal(R, ["x*y", "x^2"]).quotient("x").basis)).toEqual(["y", "x"]);
    expect(strs(ideal(R, ["x^2*y"]).saturate("x").basis)).toEqual(["y"]);
    expect(I.mul(J).isSubsetOf(I.intersect(J))).toBe(true);
    expect(I.pow(2).equals(ideal(R, ["x^2"]))).toBe(true);
    expect(I.add(ideal(R, ["x - 1"])).isWhole).toBe(true);
  });
  test("radicals", () => {
    const Q = ideal(R, ["x^2", "y^3"]);
    expect(Q.contains("x + y")).toBe(false);
    expect(Q.radicalContains("x + y")).toBe(true);
    expect(Q.exponentOf("x + y")).toBe(4);
    expect(strs(Q.radical().basis)).toEqual(["y", "x"]);
    expect(Q.isRadical()).toBe(false);
    const cert = Q.radicalContains("x + y", { certificate: true });
    expect(cert.holds).toBe(true);
    expect(cert.check).toBe(true);
    const sum = cert.cofactors.reduce((acc, c, i) => acc.add(c.mul(Q.generators[i])), R.zero);
    expect(sum.equals(R.poly("x + y").pow(cert.power))).toBe(true);
    expect(Q.radicalContains("x + 1", { certificate: true }).holds).toBe(false);
  });
  test("elimination gives the implicit equation of a curve", () => {
    const S = ring("t, x, y");
    expect(strs(ideal(S, ["x - t^2", "y - t^3"]).eliminate("t").basis)).toEqual(["x³ − y²"]);
  });
  test("prime and maximal ideals", () => {
    expect(ideal(R, ["x^2 - 2", "y - x"]).isMaximal()).toBe(true);
    expect(ideal(R, ["x^2 - 1", "y"]).isMaximal()).toBe(false);
    expect(ideal(R, ["x^2 + y^2 - 1"]).isPrime()).toBe(true);
    expect(ideal(R, ["x*y"]).isPrime()).toBe(false);
  });
  test("Chinese remainder theorem for ideals", () => {
    const A = ideal(R, ["x", "y"]);
    const B = ideal(R, ["x - 1", "y - 2"]);
    const c = chineseRemainder(A, B);
    expect(c.comaximal).toBe(true);
    expect(c.a.add(c.b).equals(R.one)).toBe(true);
    expect(c.productIsIntersection).toBe(true);
    const h = c.lift("7", "x^2 + 3");
    expect(String(A.reduce(h))).toBe("7");
    expect(String(B.reduce(h))).toBe("4");
    expect(chineseRemainder(I, J).comaximal).toBe(false);
  });
});

describe("quotient rings", () => {
  test("Q[x, y]/(x² − 2, y − x) is the field Q(√2)", () => {
    const A = ideal(R, ["x^2 - 2", "y - x"]).quotientRing();
    const x = A.gens[0];
    expect(A.dimension).toBe(2);
    expect(A.isField()).toBe(true);
    const inv = x.add(1).inverse();
    expect(inv.mul(x.add(1)).equals(A.of(1))).toBe(true);
    expect(String(x.pow(5))).toBe("4y");
    expect(strs(x.minimalPolynomial())).toEqual(["-2", "0", "1"]);
  });
  test("zero divisors, nilpotents and units", () => {
    const A = ideal(R, ["x*y", "y^2 - y"]).quotientRing();
    const y = A.gens[1];
    expect(y.isZeroDivisor()).toBe(true);
    expect(y.isUnit()).toBe(false);
    expect(A.isFinite).toBe(false);
    const B = ideal(R, ["x^2", "y"]).quotientRing();
    expect(B.gens[0].isNilpotent()).toBe(true);
    expect(B.isReduced()).toBe(false);
    expect(B.isField()).toBe(false);
    expect(B.gens[0].add(1).isUnit()).toBe(true);
    expect(() => B.gens[0].inverse()).toThrow();
  });
});

describe("solving polynomial systems", () => {
  test("rational solutions", () => {
    const s = solve(R, ["x^2 + y^2 - 5", "x*y - 2"]);
    expect(s.count).toBe(4);
    expect(strs(s.real)).toEqual(["(x = -2, y = -1)", "(x = -1, y = -2)", "(x = 1, y = 2)", "(x = 2, y = 1)"]);
  });
  test("algebraic solutions come with exact values", () => {
    const s = solve(R, ["x^2 + y^2 - 1", "y - x^2"]);
    expect(s.count).toBe(4);
    expect(s.real.length).toBe(2);
    expect(s.real[1].approx.y).toBeCloseTo((Math.sqrt(5) - 1) / 2, 12);
    expect(s.real[1].approx.x).toBeCloseTo(Math.sqrt((Math.sqrt(5) - 1) / 2), 12);
    expect(s.components.map((c) => c.degree)).toEqual([4]);
  });
  test("multiplicity, no solutions, infinitely many", () => {
    const m = solve(R, ["y - x^2", "y"]);
    expect([m.count, m.multiplicity, m.radical]).toEqual([1, 2, false]);
    const none = solve(R, ["x", "x - 1"]);
    expect(none.count).toBe(0);
    expect(solve(R, ["x*y"]).finite).toBe(false);
    const c = solve(R, ["x^2 + y^2 + 1", "x - y"], { complex: true });
    expect(c.real.length).toBe(0);
    expect(c.complex.length).toBe(2);
  });
  test("three variables", () => {
    const S = ring("x, y, z");
    const s = solve(S, ["x + y + z - 6", "x*y + y*z + z*x - 11", "x*y*z - 6"]);
    expect(s.count).toBe(6);
    expect(s.real.length).toBe(6);
    expect(s.real.every((p) => [1, 2, 3].every((v) => Object.values(p.approx).some((a) => Math.abs(a - v) < 1e-9)))).toBe(true);
  });
});

describe("limits", () => {
  test("timeLimit and maxPairs stop a large computation", () => {
    const S = ring(Array.from({ length: 7 }, (_, i) => `x${i}`), { order: "lex" });
    const g = S.gens;
    const polys = [];
    for (let k = 1; k < 7; k += 1) {
      let s = S.zero;
      for (let i = 0; i < 7; i += 1) {
        let m = S.one;
        for (let j = 0; j < k; j += 1) m = m.mul(g[(i + j) % 7]);
        s = s.add(m);
      }
      polys.push(s);
    }
    polys.push(g.reduce((m, v) => m.mul(v), S.one).sub(S.one));
    const t0 = Date.now();
    expect(() => groebner(polys, { timeLimit: 300 })).toThrow("time limit");
    expect(Date.now() - t0).toBeLessThan(5000);
    expect(() => groebner(polys, { maxPairs: 10 })).toThrow("S-polynomials");
  });
});
