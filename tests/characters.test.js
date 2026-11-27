import { test, expect } from "bun:test";
import {
  characterTable, Cyclotomic, symmetric, alternating, cyclic, dihedral, quaternion, dicyclic, PSL, GL, SL,
  metacyclic, directProduct, smallGroups, largeGroups, action,
} from "../index.js";

const rowsOf = (G) => characterTable(G).irreducibles.map((c) => c.values.map(String));

test("cyclotomic numbers: arithmetic, conjugates and printing", () => {
  const z = Cyclotomic.root(5);
  const s = z.add(z.conj());
  expect(String(s)).toBe("(−1 + √5)/2");
  expect(String(s.mul(s).add(s))).toBe("1");
  expect(s.conjugates().length).toBe(2);
  expect(s.conductor).toBe(5);
  expect(String(Cyclotomic.root(4))).toBe("i");
  expect(String(Cyclotomic.root(3))).toBe("(−1 + i√3)/2");
  expect(String(Cyclotomic.root(8).add(Cyclotomic.root(8, 3)))).toBe("i√2");
  expect(String(Cyclotomic.root(8).add(Cyclotomic.root(8, 7)))).toBe("√2");
  expect(String(Cyclotomic.root(7, 6))).toBe("ζ₇⁶");
  expect(String(Cyclotomic.root(12, 4))).toBe("(−1 + i√3)/2");
  expect(Cyclotomic.root(12, 4).conductor).toBe(3);
  expect(Cyclotomic.root(6).lift(12).equals(Cyclotomic.root(12, 2))).toBe(true);
  expect(Cyclotomic.root(3).add(Cyclotomic.root(4)).n).toBe(12);
  expect(Cyclotomic.int(7).over(2).toFraction()).toEqual([7n, 2n]);
  expect(Cyclotomic.root(5).add(Cyclotomic.root(5, 4)).toTeX()).toBe("\\tfrac{-1 + \\sqrt{5}}{2}");
  const a = Cyclotomic.root(9).approx();
  expect(a.re).toBeCloseTo(Math.cos((2 * Math.PI) / 9), 12);
  expect(a.im).toBeCloseTo(Math.sin((2 * Math.PI) / 9), 12);
});

test("the classical tables", () => {
  expect(rowsOf(symmetric(3))).toEqual([["1", "1", "1"], ["1", "−1", "1"], ["2", "0", "−1"]]);
  expect(rowsOf(symmetric(4))).toEqual([
    ["1", "1", "1", "1", "1"],
    ["1", "1", "−1", "1", "−1"],
    ["2", "2", "0", "−1", "0"],
    ["3", "−1", "1", "0", "−1"],
    ["3", "−1", "−1", "0", "1"],
  ]);
  expect(rowsOf(alternating(5))).toEqual([
    ["1", "1", "1", "1", "1"],
    ["3", "−1", "0", "(1 + √5)/2", "(1 − √5)/2"],
    ["3", "−1", "0", "(1 − √5)/2", "(1 + √5)/2"],
    ["4", "0", "1", "−1", "−1"],
    ["5", "1", "−1", "0", "0"],
  ]);
  expect(rowsOf(cyclic(4))[2]).toEqual(["1", "−1", "i", "−i"]);
  const T = characterTable(PSL(2, 7));
  expect(T.degrees).toEqual([1, 3, 3, 6, 7, 8]);
  expect(String(T.irreducibles[1].values[4])).toBe("(−1 + i√7)/2");
});

test("D4 and Q8 have the same table but different indicators", () => {
  const D = characterTable(dihedral(4));
  const Q = characterTable(quaternion());
  expect(D.irreducibles.map((c) => c.values.map(String))).toEqual(Q.irreducibles.map((c) => c.values.map(String)));
  expect(D.irreducibles.map((c) => c.indicator())).toEqual([1, 1, 1, 1, 1]);
  expect(Q.irreducibles.map((c) => c.indicator())).toEqual([1, 1, 1, 1, -1]);
  expect(Q.irreducibles[4].type).toBe("quaternionic");
  expect(characterTable(cyclic(3)).irreducibles[1].type).toBe("complex");
});

test("orthogonality, degrees and counts for a range of groups", () => {
  const groups = [cyclic(12), dihedral(6), dihedral(7), dicyclic(3), dicyclic(5), symmetric(5), alternating(6), GL(2, 3), SL(2, 3), SL(2, 5),
    metacyclic(7, 3, 2), directProduct(symmetric(3), cyclic(4)), ...smallGroups(16), ...smallGroups(18), ...smallGroups(20)];
  for (const G of groups) {
    const T = characterTable(G);
    expect(T.verify()).toEqual({ rows: true, columns: true, degrees: true, ok: true });
    expect(T.irreducibles.length).toBe(G.conjugacyClasses().length);
    expect(T.degrees.reduce((s, d) => s + d * d, 0)).toBe(G.order);
    expect(T.abelianizationOrder).toBe(G.order / G.derivedSubgroup().order);
    expect(T.centreOrder).toBe(G.centre().order);
    expect(T.normalSubgroups().map((N) => N.order)).toEqual(G.normalSubgroups().map((N) => N.order).sort((a, b) => a - b));
    // real characters against real classes
    expect(T.irreducibles.filter((c) => c.values.every((v) => v.isReal)).length).toBe(T.realClasses.length);
    // Σ indicator·degree = number of square roots of the identity
    const inv = G.elements.filter((x) => G.elementOrder(x) <= 2).length;
    expect(T.irreducibles.reduce((s, c) => s + c.indicator() * c.degree, 0)).toBe(inv);
  }
});

test("larger tables: S6, S7 and M11", () => {
  expect(characterTable(symmetric(6)).degrees).toEqual([1, 1, 5, 5, 5, 5, 9, 9, 10, 10, 16]);
  expect(characterTable(symmetric(7)).degrees).toEqual([1, 1, 6, 6, 14, 14, 14, 14, 15, 15, 20, 21, 21, 35, 35]);
  const T = characterTable(largeGroups.mathieu(11).toGroup());
  expect(T.degrees).toEqual([1, 10, 10, 10, 11, 16, 16, 44, 45, 55]);
  expect(T.isSimple).toBe(true);
  expect(T.irreducibles[5].values.map(String).slice(-2)).toEqual(["(−1 + i√11)/2", "(−1 − i√11)/2"]);
});

test("products, squares and decomposition", () => {
  const T = characterTable(symmetric(4));
  const [, sign, two, std, other] = T.irreducibles;
  expect(std.mul(sign).equals(other)).toBe(true);
  expect(std.mul(std).decomposition()).toBe("χ₁ + χ₃ + χ₄ + χ₅");
  expect(std.symmetricSquare().decomposition()).toBe("χ₁ + χ₃ + χ₄");
  expect(std.alternatingSquare().decomposition()).toBe("χ₅");
  expect(T.regular.multiplicities()).toEqual(T.degrees);
  expect(two.kernel().order).toBe(4);
  expect(std.isFaithful).toBe(true);
  expect(std.isIrreducible).toBe(true);
  expect(std.add(two).isIrreducible).toBe(false);
  expect(std.add(two).isCharacter).toBe(true);
  expect(std.sub(two).isCharacter).toBe(false);
  expect(std.sub(two).isVirtualCharacter).toBe(true);
  expect(std.inner(std)).toBe(1);
  expect(std.inner(two)).toBe(0);
  expect(std.adams(2).decomposition()).toBe("χ₁ + χ₃ + χ₄ − χ₅");
  expect(T.classFunction([4, 0, 2, 1, 0]).decomposition()).toBe("χ₁ + χ₄");
});

test("restriction, induction and Frobenius reciprocity", () => {
  const G = symmetric(4);
  const T = characterTable(G);
  for (const H of [G.sylowSubgroup(2), G.sylowSubgroup(3), G.derivedSubgroup()]) {
    const TH = characterTable(H);
    for (const chi of T.irreducibles) for (const phi of TH.irreducibles) expect(phi.induce(T).inner(chi)).toBe(chi.restrict(TH).inner(phi));
    expect(TH.trivial.induce(T).degree).toBe(G.order / H.order);
  }
  // the action on four points: 1 + standard, so it is doubly transitive
  const X = action(G, [0, 1, 2, 3], (g, x) => g[x]);
  const pi = T.permutationCharacter(X);
  expect(pi.values.map(String)).toEqual(["4", "0", "2", "1", "0"]);
  expect(pi.norm).toBe(2);
  expect(pi.decomposition()).toBe("χ₁ + χ₄");
  expect(T.permutationCharacter(G.sylowSubgroup(2)).decomposition()).toBe("χ₁ + χ₃");
});

test("printing", () => {
  const T = characterTable(symmetric(3));
  expect(T.toString()).toBe(["      1a  2a  3a", "size   1   3   2", "χ₁     1   1   1", "χ₂     1  −1   1", "χ₃     2   0  −1"].join("\n"));
  expect(T.toTeX()).toContain("\\chi_{3} & 2 & 0 & -1");
  expect(T.describe()).toBe("S3 of order 6: 3 classes, degrees 1, 1, 2; values in Q(ζ₆), found mod 7");
  expect(T.classes.map((c) => c.name)).toEqual(["1a", "2a", "3a"]);
  expect(T.powerMap(2)).toEqual([0, 0, 2]);
  const G = symmetric(3);
  expect(characterTable(G) === characterTable(G)).toBe(true);
});
