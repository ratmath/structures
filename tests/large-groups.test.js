import { describe, expect, test } from "bun:test";
import { schreierSims, largeGroups, symmetric, dihedral, GL, PermGroup } from "../index.js";

const factorial = (n) => { let f = 1n; for (let i = 2n; i <= BigInt(n); i += 1n) f *= i; return f; };

describe("Schreier–Sims", () => {
  test("orders agree with groups that can be listed", () => {
    for (const n of [1, 2, 3, 4, 5, 6, 7]) {
      expect(largeGroups.symmetric(n).order).toBe(factorial(n));
      expect(largeGroups.alternating(n).order).toBe(n < 2 ? 1n : factorial(n) / 2n);
    }
    const D = schreierSims(["(1 2 3 4 5 6)", "(2 6)(3 5)"]);
    expect(D.order).toBe(12n);
    expect(D.toGroup().order).toBe(dihedral(6).order);
    // GL(2, 3) acting on the 8 non-zero vectors
    const G = GL(2, 3);
    const vectors = [];
    for (let a = 0; a < 3; a += 1) for (let b = 0; b < 3; b += 1) if (a || b) vectors.push([a, b]);
    const at = (v) => vectors.findIndex((w) => w[0] === v[0] && w[1] === v[1]);
    const perm = (m) => vectors.map(([x, y]) => at([Number((BigInt(m[0][0]) * BigInt(x) + BigInt(m[0][1]) * BigInt(y)) % 3n), Number((BigInt(m[1][0]) * BigInt(x) + BigInt(m[1][1]) * BigInt(y)) % 3n)]));
    expect(schreierSims(G.generators.map(perm), { degree: 8 }).order).toBe(48n);
  });
  test("the order is the product of the basic orbits", () => {
    const S = largeGroups.symmetric(50);
    expect(S.order).toBe(factorial(50));
    expect(S.orbitLengths.reduce((a, b) => a * BigInt(b), 1n)).toBe(S.order);
    expect(S.isSymmetric).toBe(true);
    expect(largeGroups.alternating(30).isAlternating).toBe(true);
    expect(largeGroups.alternating(30).isSymmetric).toBe(false);
  });
  test("the Mathieu groups have their known orders", () => {
    expect(largeGroups.mathieu(11).order).toBe(7920n);
    expect(largeGroups.mathieu(12).order).toBe(95040n);
    expect(largeGroups.mathieu(22).order).toBe(443520n);
    expect(largeGroups.mathieu(23).order).toBe(10200960n);
    expect(largeGroups.mathieu(24).order).toBe(244823040n);
    const M12 = largeGroups.mathieu(12);
    expect(M12.isPrimitive).toBe(true);
    expect(M12.isPerfect).toBe(true);
    // sharply 5-transitive: fixing five points leaves only the identity
    expect(M12.stabilizer([1, 2, 3, 4, 5]).order).toBe(1n);
    expect(M12.stabilizer([1]).order).toBe(7920n);
    expect(() => largeGroups.mathieu(13)).toThrow();
  });
  test("Rubik's cube", () => {
    const R = largeGroups.rubik();
    expect(R.order).toBe(43252003274489856000n);
    expect(R.degree).toBe(48);
    expect(R.orbits().map((o) => o.length)).toEqual([24, 24]);
    // a single swap of two stickers cannot be reached
    expect(R.contains("(1 2)")).toBe(false);
    expect(R.sift("(1 2)").member).toBe(false);
    const g = R.evaluate(["R", "U", "R⁻¹", "U⁻¹"]);
    expect(R.contains(g)).toBe(true);
    const w = R.word(g);
    expect(w.check).toBe(true);
    expect(R.evaluate(w.word)).toEqual(g);
    // the commutator subgroup has index 2
    expect(R.index(R.derivedSubgroup())).toBe(2n);
    expect(R.isSolvable).toBe(false);
    expect(largeGroups.pocketCube().order).toBe(3674160n);
  });
  test("membership, words and random elements", () => {
    const S = largeGroups.symmetric(8);
    const g = S.random({ seed: 7 });
    expect(S.random({ seed: 7 })).toEqual(g);
    expect(S.contains(g)).toBe(true);
    const w = S.word(g);
    expect(w.check).toBe(true);
    expect(S.sift(g).member).toBe(true);
    const A = largeGroups.alternating(8);
    expect(A.contains("(1 2)")).toBe(false);
    expect(A.contains("(1 2 3)")).toBe(true);
    expect(S.isNormal(A)).toBe(true);
    expect(S.index(A)).toBe(2n);
    expect([...largeGroups.symmetric(4).elements()].length).toBe(24);
    expect(new Set([...largeGroups.symmetric(4).elements()].map(String)).size).toBe(24);
  });
  test("series, closures and blocks", () => {
    const S4 = largeGroups.symmetric(4);
    expect(S4.derivedSeries().map((G) => G.order)).toEqual([24n, 12n, 4n, 1n]);
    expect(S4.isSolvable).toBe(true);
    expect(largeGroups.symmetric(5).isSolvable).toBe(false);
    expect(S4.normalClosure(["(1 2)(3 4)"]).order).toBe(4n);
    const D = schreierSims(["(1 2 3 4 5 6)", "(2 6)(3 5)"]);
    expect(D.blockSystems().map((b) => b.length).sort()).toEqual([2, 3]);
    expect(D.isPrimitive).toBe(false);
    expect(D.minimalBlocks(1, 4)).toEqual([[1, 4], [2, 5], [3, 6]]);
    expect(largeGroups.symmetric(6).isPrimitive).toBe(true);
    expect(D.isAbelian).toBe(false);
    expect(schreierSims(["(1 2 3)", "(4 5)"]).isAbelian).toBe(true);
    expect(schreierSims(["(1 2 3)", "(4 5)"]).isTransitive).toBe(false);
  });
  test("a prescribed base, and subgroups", () => {
    const S = schreierSims(["(1 2)", "(1 2 3 4 5)"], { base: [5, 4] });
    expect(S.base.slice(0, 2)).toEqual([5, 4]);
    expect(S.order).toBe(120n);
    const H = S.subgroup(["(1 2 3)", "(1 2)"]);
    expect(H.order).toBe(6n);
    expect(S.isSubgroup(H)).toBe(true);
    expect(() => largeGroups.alternating(5).subgroup(["(1 2)"])).toThrow();
    expect(schreierSims(new PermGroup([[1, 2, 0]], 3)).order).toBe(3n);
    expect(S.equals(largeGroups.symmetric(5))).toBe(true);
    expect(symmetric(5).order).toBe(120);
  });
});
