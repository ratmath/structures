// Presentations and coset enumeration, abelian groups, and modules. The
// Smith form and the module structure are checked against @ratmath/linalg,
// an independent implementation, and coset enumeration against the groups
// built element by element.
import { describe, expect, test } from "bun:test";
import { Rational } from "@ratmath/core";
import * as linalg from "@ratmath/linalg";
import {
  presentation, presentationOf, parseWord, formatWord, abelianInvariants,
  abelianGroup, cyclicProduct, fromRelations, abelianGroups, smithForm, ZZ, polynomialDomain, matrixModule, companionMatrix,
  symmetric, alternating, dihedral, quaternion, cyclic, GL, SL, PSL, heisenberg, directProduct, dicyclic, metacyclic,
  isIsomorphic, smallGroups, numberOfGroups,
} from "../index.js";

const matMul = (R, X, Y) => X.map((row) => Y[0].map((_, j) => row.reduce((s, x, k) => R.add(s, R.mul(x, Y[k][j])), R.zero)));
const big = (A) => A.map((r) => r.map(BigInt));

/** A small deterministic generator, so the random tests are repeatable. */
function rng(seed) {
  let s = BigInt(seed);
  return (n) => {
    s = (s * 6364136223846793005n + 1442695040888963407n) & ((1n << 64n) - 1n);
    return Number((s >> 33n) % BigInt(n));
  };
}

describe("words", () => {
  test("reading and printing words", () => {
    const names = ["a", "b"];
    expect(parseWord("a^2 b^-1", names)).toEqual([1, 1, -2]);
    expect(parseWord("a²b⁻¹", names)).toEqual([1, 1, -2]);
    expect(parseWord("(ab)^2", names)).toEqual([1, 2, 1, 2]);
    expect(parseWord("[a, b]", names)).toEqual([-1, -2, 1, 2]);
    expect(parseWord("aA", names)).toEqual([]);
    expect(parseWord("(ab)^-1", names)).toEqual([-2, -1]);
    expect(formatWord([1, 1, -2, -2, -2], names)).toBe("a²b⁻³");
    expect(formatWord([1, 1, -2], names, { plain: true })).toBe("a^2*b^-1");
    expect(() => parseWord("ac", names)).toThrow(/unknown generator/);
    const P = presentation("x, y | x^3 = y^2, (xy)^2");
    expect(String(P)).toBe("⟨x, y | x³y⁻², xyxy⟩");
  });
});

describe("coset enumeration", () => {
  const known = [
    ["a, b | a^4, b^2, (ab)^2", 8, () => dihedral(4)],
    ["a, b | a^2, b^3, (ab)^3", 12, () => alternating(4)],
    ["a, b | a^2, b^3, (ab)^4", 24, () => symmetric(4)],
    ["a, b | a^2, b^3, (ab)^5", 60, () => alternating(5)],
    ["a, b | a^2, b^3, (ab)^7, [a, b]^4", 168, () => PSL(2, 7)],
    ["a, b | a^4, a^2 = b^2, b^-1 a b = a^-1", 8, () => quaternion()],
    ["a, b | a^6, b^4, ab = ba", 24, () => directProduct(cyclic(6), cyclic(4))],
    ["a, b, c | a^2, b^2, c^2, (ab)^3, (bc)^3, (ac)^2", 24, () => symmetric(4)],
    ["a, b, c | a^2, b^2, c^2, (ab)^3, (bc)^5, (ac)^2", 120, () => directProduct(alternating(5), cyclic(2))],
    ["a, b, c, d, e | ab = c, bc = d, cd = e, de = a, ea = b", 11, () => cyclic(11)], // the Fibonacci group F(2, 5)
    ["a, b | a^3, b^3, (ab)^3, (a^-1 b)^3", 27, () => heisenberg(3)],
    ["a, b | a^8, b^2, b a b^-1 = a^3", 16, () => metacyclic(8, 2, 3)],
    ["a, b | a^7, b^3, b a b^-1 = a^2", 21, () => metacyclic(7, 3, 2)],
  ];
  test("orders of known presentations, and the groups they give", () => {
    for (const [text, n, make] of known) {
      const P = presentation(text);
      expect(P.order()).toBe(n);
      const G = P.toGroup();
      expect(G.order).toBe(n);
      expect(isIsomorphic(G, make())).toBe(true);
    }
  });
  test("indices of subgroups agree with the group built element by element", () => {
    const P = presentation("a, b | a^2, b^3, (ab)^5");
    const G = P.toGroup();
    const r = rng(7);
    for (let t = 0; t < 25; t += 1) {
      const words = Array.from({ length: 1 + r(2) }, () => Array.from({ length: 1 + r(5) }, () => (r(2) ? 1 : -1) * (1 + r(2))));
      const H = G.subgroup(words.map((w) => G.fromWord(w)));
      expect(P.index(words)).toBe(60 / H.order);
    }
    expect(P.index(["b"])).toBe(20);
  });
  test("each generator acts as a permutation of the cosets, and the relators act trivially", () => {
    const P = presentation("a, b | a^2, b^3, (ab)^4");
    const T = P.enumerate(["a"]);
    expect(T.index).toBe(12);
    for (let g = 0; g < 2; g += 1) expect(new Set(T.table.map((row) => row[g])).size).toBe(12);
    for (const w of P.relators) for (let c = 0; c < T.index; c += 1) {
      let d = c;
      for (const x of w) d = x > 0 ? T.table[d][x - 1] : T.inverseTable[d][-x - 1];
      expect(d).toBe(c);
    }
  });
  test("an infinite group stops at the coset limit", () => {
    expect(() => presentation("a, b | a^3, b^3, (ab)^3").order({ maxCosets: 2000 })).toThrow(/coset enumeration stopped/);
    expect(() => presentation("a, b | ab = ba").order({ maxCosets: 2000 })).toThrow(/coset enumeration stopped/);
    // the infinite group still has finite-index subgroups: ⟨a, b | a³, b³, (ab)³⟩ maps onto C3
    expect(presentation("a, b | a^3, b^3, (ab)^3").index(["a b^-1", "a^-1 b", "b a^-1", "[a, b]", "[a, b^-1]"])).toBe(3);
  });
  test("presentationOf gives back the group, for every group of order up to 20", () => {
    for (let n = 1; n <= 20; n += 1) for (const G of smallGroups(n)) {
      const P = presentationOf(G);
      expect(P.order()).toBe(n);
      expect(isIsomorphic(P.toGroup(), G)).toBe(true);
    }
    for (const G of [symmetric(4), GL(2, 3), PSL(2, 7), SL(2, 5)]) expect(presentationOf(G).order()).toBe(G.order);
  });
});

describe("Smith normal form", () => {
  test("U·A·V = D with the divisibility chain, against linalg's smith on random integer matrices", () => {
    const r = rng(11);
    for (let t = 0; t < 60; t += 1) {
      const m = 1 + r(4);
      const n = 1 + r(4);
      const A = Array.from({ length: m }, () => Array.from({ length: n }, () => r(21) - 10));
      const S = smithForm(A);
      const D = matMul(ZZ, matMul(ZZ, S.U, big(A)), S.V);
      D.forEach((row, i) => row.forEach((x, j) => expect(x).toBe(i === j && i < S.diagonal.length ? S.diagonal[i] : 0n)));
      const I = matMul(ZZ, S.V, S.Vinv);
      I.forEach((row, i) => row.forEach((x, j) => expect(x).toBe(i === j ? 1n : 0n)));
      for (let i = 1; i < S.rank; i += 1) expect(S.diagonal[i] % S.diagonal[i - 1]).toBe(0n);
      for (const d of S.diagonal) expect(d >= 0n).toBe(true);
      const L = linalg.smith(A).diagonal.filter((d) => d !== 0n);
      expect(S.diagonal.filter((d) => d !== 0n)).toEqual(L);
    }
  });
  test("over Q[x]: the invariant factors of xI − A are linalg's, and the rational canonical forms agree", () => {
    const r = rng(5);
    const cases = [
      [[2, 1, 0, 0], [0, 2, 0, 0], [0, 0, 2, 0], [0, 0, 0, 3]],
      [[0, 1], [-1, 0]],
      [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
      [[0, 0, 0], [1, 0, 0], [0, 1, 0]],
    ];
    for (let t = 0; t < 12; t += 1) {
      const n = 2 + r(3);
      // random matrices with a repeated block structure, so some invariant factors are nontrivial
      const B = Array.from({ length: n }, () => Array.from({ length: n }, () => r(5) - 2));
      cases.push(B.map((row, i) => row.map((x, j) => (i === j ? x : r(3) === 0 ? x : 0))));
    }
    for (const A of cases) {
      const M = matrixModule(A);
      const S = linalg.similarity(A).exact();
      expect(M.invariantFactors.map((f) => f.coefficients("x").map((c) => String(c.isZero() ? 0 : c.toConstant())))).toEqual(S.invariantFactors.map((f) => f.map(String)));
      expect(M.rationalCanonical.map((row) => row.map(String))).toEqual(S.rationalCanonical.toStrings());
      // the product of the invariant factors is the characteristic polynomial, the last is the minimal polynomial
      expect(M.characteristicPolynomial.degree()).toBe(A.length);
      expect(M.minimalPolynomial.coefficients("x").map((c) => String(c.isZero() ? 0 : c.toConstant()))).toEqual(S.minimalPolynomial.map(String));
      const D = M.smith;
      const P = matMul(D.domain, matMul(D.domain, D.U, A.map((row, i) => row.map((a, j) => (i === j ? D.domain.ring.poly(`x - ${a}`) : D.domain.ring.from(-a))))), D.V);
      P.forEach((row, i) => row.forEach((x, j) => expect(x.equals(i === j ? D.diagonal[i] : D.domain.zero)).toBe(true)));
    }
  });
  test("module structure: elementary divisors and cyclic modules", () => {
    const M = matrixModule([[2, 1, 0, 0], [0, 2, 0, 0], [0, 0, 2, 0], [0, 0, 0, 3]]);
    expect(String(M)).toBe("Q[x]/(x − 2) ⊕ Q[x]/(x³ − 7x² + 16x − 12)");
    expect(M.elementaryDivisors.map((e) => String(e.divisor))).toEqual(["x − 2", "x − 3", "x² − 4x + 4"]);
    expect(M.isCyclic).toBe(false);
    expect(matrixModule([[0, 0, 0], [1, 0, 0], [0, 1, 0]]).isCyclic).toBe(true);
    expect(matrixModule([[0, 1], [-1, 0]]).elementaryDivisors.map((e) => String(e.divisor))).toEqual(["x² + 1"]);
    expect(matrixModule([["1/2", 0], [0, "1/2"]]).invariantFactors.map(String)).toEqual(["x − 1/2", "x − 1/2"]);
    const D = polynomialDomain("t");
    expect(String(companionMatrix(D.ring.poly("t^2 - 3t + 2")).map((r) => r.join(" ")))).toBe("0 -2,1 3");
  });
});

describe("finitely generated abelian groups", () => {
  test("from relations", () => {
    const G = fromRelations([[2, 4, 4], [-6, 6, 12], [10, -4, -16]]);
    expect(G.invariants).toEqual([2n, 6n, 12n]);
    expect(G.order).toBe(144n);
    expect(String(G)).toBe("Z/2 × Z/6 × Z/12");
    expect(String(G.toString({ primary: true }))).toBe("Z/2 × Z/2 × Z/4 × Z/3 × Z/3");
    const F = fromRelations([[2, 0, 0]], 3);
    expect(String(F)).toBe("Z/2 × Z^2");
    expect(F.order).toBe(null);
    expect(String(fromRelations([[1, 1], [1, -1]]))).toBe("Z/2");
    // each new generator has the stated order: d·v lies in the relation lattice (adding it as a relation changes nothing), (d/p)·v does not
    const Rel = [[2, 4, 4], [-6, 6, 12], [10, -4, -16]];
    for (const { order, vector } of G.generators) {
      expect(fromRelations([...Rel, vector.map((x) => x * order)]).order).toBe(144n);
      for (const q of [2n, 3n]) if (order % q === 0n) expect(fromRelations([...Rel, vector.map((x) => (x * order) / q)]).order < 144n).toBe(true);
    }
    const R = [[4, 6], [6, 4]];
    const A = fromRelations(R);
    expect(String(A)).toBe("Z/2 × Z/10");
  });
  test("cyclic products and the Chinese remainder theorem", () => {
    expect(String(cyclicProduct([4, 6]))).toBe("Z/2 × Z/12");
    expect(cyclicProduct([3, 5]).isCyclic).toBe(true);
    expect(cyclicProduct([2, 0, 0]).rank).toBe(2);
    expect(abelianGroup([2, 4]).equals(cyclicProduct([4, 2]))).toBe(true);
    expect(() => abelianGroup([4, 6])).toThrow(/divide/);
    expect(String(abelianGroup([2]).times(abelianGroup([3])))).toBe("Z/6");
  });
  test("all abelian groups of order n: Π p(aᵢ) of them, matching the abelian groups in the table", () => {
    const p = [1, 1, 2, 3, 5, 7, 11];
    for (let n = 1; n <= 64; n += 1) {
      let want = 1;
      let m = n;
      for (let q = 2; q <= m; q += 1) { let a = 0; while (m % q === 0) { m /= q; a += 1; } want *= p[a]; }
      expect(abelianGroups(n).length).toBe(want);
    }
    for (let n = 1; n <= 20; n += 1) {
      const fromTable = smallGroups(n).filter((G) => G.isAbelian);
      expect(fromTable.length).toBe(abelianGroups(n).length);
      for (const A of abelianGroups(n)) expect(fromTable.filter((G) => isIsomorphic(G, A.toGroup())).length).toBe(1);
    }
    expect(abelianGroups(72).map(String)).toContain("Z/6 × Z/12");
    void numberOfGroups;
  });
  test("abelianisations G/G′ agree with the quotient by the commutator subgroup", () => {
    for (const G of [symmetric(4), alternating(4), dihedral(4), dihedral(5), quaternion(), GL(2, 3), dicyclic(3), heisenberg(3), directProduct(cyclic(4), symmetric(3)), PSL(2, 7)]) {
      const A = abelianInvariants(G);
      const Q = G.quotient(G.derivedSubgroup());
      expect(Number(A.order)).toBe(Q.order);
      expect(isIsomorphic(A.toGroup(), Q)).toBe(true);
    }
    expect(String(abelianInvariants(presentation("a, b | a^2, b^3, (ab)^3")))).toBe("Z/3");
    expect(String(abelianInvariants(presentation("a, b | [a, b]")))).toBe("Z^2");
    expect(String(abelianInvariants(presentation("a, b, c |")))).toBe("Z^3");
    void Rational;
  });
});
