// Finite groups: every algorithm is checked against brute force on small
// groups (working only with the elements and the multiplication) and
// against known facts.
import { describe, expect, test } from "bun:test";
import {
  cyclic, dihedral, dicyclic, quaternion, kleinFour, symmetric, alternating, permutationGroup, fromPermGroup,
  matrixGroup, GL, SL, PSL, heisenberg, directProduct, semidirectProduct, metacyclic, Group, PermGroup,
} from "../index.js";

// ---- brute force, from the multiplication alone ---------------------------------------------------

const keys = (G, xs) => new Set(xs.map((x) => G.key(x)));
const sameSet = (G, a, b) => {
  const A = keys(G, a);
  const B = keys(G, b);
  return A.size === B.size && [...A].every((k) => B.has(k));
};
const bruteOrder = (G, x) => {
  let k = 1;
  for (let y = x; G.key(y) !== G.key(G.identity); y = G.mul(y, x)) k += 1;
  return k;
};
const bruteCentre = (G) => G.elements.filter((z) => G.elements.every((g) => G.key(G.mul(z, g)) === G.key(G.mul(g, z))));
const bruteClass = (G, x) => G.elements.map((g) => G.mul(G.mul(g, x), G.inv(g)));
const bruteClose = (G, gens) => {
  const out = new Map([[G.key(G.identity), G.identity]]);
  let grew = true;
  while (grew) {
    grew = false;
    for (const a of [...out.values()]) for (const b of gens) {
      const c = G.mul(a, b);
      if (!out.has(G.key(c))) { out.set(G.key(c), c); grew = true; }
    }
  }
  return [...out.values()];
};
const bruteIsNormal = (G, H) => {
  const hk = keys(G, H.elements);
  return G.elements.every((g) => H.elements.every((h) => hk.has(G.key(G.mul(G.mul(g, h), G.inv(g))))));
};
const isClosed = (G, xs) => {
  const k = keys(G, xs);
  return xs.every((a) => xs.every((b) => k.has(G.key(G.mul(a, b)))));
};
const divisors = (n) => [...Array(n).keys()].map((i) => i + 1).filter((d) => n % d === 0);
const phi = (n) => [...Array(n).keys()].filter((i) => gcdN(i, n) === 1).length;
function gcdN(a, b) { while (b) [a, b] = [b, a % b]; return a; }
const sigma = (n) => divisors(n).reduce((s, d) => s + d, 0);

const SMALL = () => [
  cyclic(1), cyclic(2), cyclic(7), cyclic(12), kleinFour(), dihedral(3), dihedral(4), dihedral(5), dihedral(6), quaternion(), dicyclic(3),
  symmetric(3), symmetric(4), alternating(4), GL(2, 2), SL(2, 3), heisenberg(3), directProduct(cyclic(2), dihedral(4)), metacyclic(5, 4, 2),
];
const MEDIUM = () => [alternating(5), GL(2, 3), symmetric(5)];

describe("the group axioms and element orders", () => {
  test("every named group is closed, has inverses and is associative", () => {
    for (const G of [...SMALL(), GL(2, 3)]) {
      const E = G.elements;
      expect(isClosed(G, E)).toBe(true);
      for (const x of E) expect(G.key(G.mul(x, G.inv(x)))).toBe(G.key(G.identity));
      for (const a of E.slice(0, 12)) for (const b of E.slice(0, 12)) for (const c of E.slice(0, 12)) expect(G.key(G.mul(G.mul(a, b), c))).toBe(G.key(G.mul(a, G.mul(b, c))));
    }
  });
  test("orders of the families", () => {
    for (let n = 1; n <= 10; n += 1) {
      expect(cyclic(n).order).toBe(n);
      expect(dihedral(n).order).toBe(2 * n);
      expect(dicyclic(n).order).toBe(4 * n);
    }
    expect([1, 2, 3, 4, 5].map((n) => symmetric(n).order)).toEqual([1, 2, 6, 24, 120]);
    expect([1, 2, 3, 4, 5].map((n) => alternating(n).order)).toEqual([1, 1, 3, 12, 60]);
    // |GL(n,p)| = (pⁿ − 1)(pⁿ − p)…(pⁿ − pⁿ⁻¹); |SL| = |GL|/(p − 1)
    for (const [n, p] of [[1, 5], [2, 2], [2, 3], [2, 5], [3, 2]]) {
      let o = 1;
      for (let i = 0; i < n; i += 1) o *= p ** n - p ** i;
      expect(GL(n, p).order).toBe(o);
      expect(SL(n, p).order).toBe(o / (p - 1));
    }
    expect(PSL(2, 5).order).toBe(60);
    expect(PSL(2, 7).order).toBe(168);
    expect(heisenberg(5).order).toBe(125);
  });
  test("element orders against brute force, and Lagrange for elements", () => {
    for (const G of SMALL()) {
      for (const x of G.elements) {
        expect(G.elementOrder(x)).toBe(bruteOrder(G, x));
        expect(G.order % G.elementOrder(x)).toBe(0);
      }
      expect([...G.orderStatistics().values()].reduce((a, b) => a + b, 0)).toBe(G.order);
    }
    expect(quaternion().orderStatistics()).toEqual(new Map([[1, 1], [2, 1], [4, 6]]));
    expect(symmetric(4).exponent).toBe(12);
  });
  test("elements can be named by their printed labels", () => {
    const D = dihedral(4);
    expect(D.mul("r", "s")).toEqual([1, 1]);
    expect(D.format(D.mul("s", "r"))).toBe("r³s");
    expect(D.format(D.pow("r", -1))).toBe("r³");
    const S = symmetric(3);
    expect(S.format(S.mul("(1 2)", "(2 3)"))).toBe("(1 2 3)");
    const Q = quaternion();
    expect(Q.format(Q.mul("i", "j"))).toBe("k");
    expect(Q.format(Q.mul("j", "i"))).toBe("−k");
    expect(Q.format(Q.mul("-1", "-1"))).toBe("1");
    expect(() => D.element("t")).toThrow(/not an element/);
  });
  test("fromTable checks the axioms", () => {
    const G = Group.fromTable([[0, 1, 2], [1, 2, 0], [2, 0, 1]]);
    expect(G.isCyclic).toBe(true);
    expect(() => Group.fromTable([[0, 1, 2], [1, 0, 0], [2, 1, 1]])).toThrow();
    // a Latin square with an identity that is not associative (a loop, not a group)
    const loop = [[0, 1, 2, 3, 4], [1, 0, 3, 4, 2], [2, 4, 0, 1, 3], [3, 2, 4, 0, 1], [4, 3, 1, 2, 0]];
    expect(() => Group.fromTable(loop)).toThrow(/associative/);
  });
  test("Cayley table rows are permutations (the Latin square property)", () => {
    const T = dihedral(5).cayleyTable();
    for (const r of T.table) expect(new Set(r).size).toBe(10);
    for (let j = 0; j < 10; j += 1) expect(new Set(T.table.map((r) => r[j])).size).toBe(10);
    expect(String(cyclic(3).cayleyTable()).split("\n")[2].trim()).toBe("e  │ e  a  a²");
  });
});

describe("centre, classes, commutators", () => {
  test("centre and centralisers against brute force", () => {
    for (const G of [...SMALL(), ...MEDIUM()]) {
      expect(sameSet(G, G.centre().elements, bruteCentre(G))).toBe(true);
      for (const x of G.elements.slice(0, 6)) {
        const C = G.centraliser(x);
        expect(sameSet(G, C.elements, G.elements.filter((g) => G.key(G.mul(g, x)) === G.key(G.mul(x, g))))).toBe(true);
      }
    }
  });
  test("conjugacy classes against brute force; the class equation", () => {
    for (const G of [...SMALL(), ...MEDIUM()]) {
      const cls = G.conjugacyClasses();
      expect(cls.reduce((s, c) => s + c.size, 0)).toBe(G.order);
      for (const c of cls) {
        expect(sameSet(G, c.elements, bruteClass(G, c.representative))).toBe(true);
        expect(c.size * G.centraliser(c.representative).order).toBe(G.order); // orbit–stabiliser
      }
      const eq = G.classEquation();
      expect(eq.centre + eq.terms.reduce((a, b) => a + b, 0)).toBe(G.order);
      expect(eq.centre).toBe(G.centre().order);
    }
  });
  test("known numbers of classes", () => {
    // Sₙ: the partitions of n; A₅: 5; GL(2, p): p² − 1; Q₈ and D₄: 5
    expect([1, 2, 3, 4, 5].map((n) => symmetric(n).conjugacyClasses().length)).toEqual([1, 2, 3, 5, 7]);
    expect(alternating(5).conjugacyClasses().length).toBe(5);
    for (const p of [2, 3, 5]) expect(GL(2, p).conjugacyClasses().length).toBe(p * p - 1);
    expect(String(alternating(5).classEquation())).toBe("60 = 1 + 12 + 12 + 15 + 20");
    // a p-group has a nontrivial centre
    for (const G of [dihedral(4), quaternion(), heisenberg(3), dihedral(8)]) expect(G.centre().order).toBeGreaterThan(1);
  });
  test("commutator subgroup against brute force", () => {
    for (const G of [...SMALL(), alternating(5), GL(2, 3)]) {
      const comms = [];
      for (const a of G.elements) for (const b of G.elements) comms.push(G.commutator(a, b));
      expect(sameSet(G, G.derivedSubgroup().elements, bruteClose(G, comms))).toBe(true);
      // G/G′ is abelian
      expect(G.quotient(G.derivedSubgroup()).isAbelian).toBe(true);
    }
    for (const n of [2, 3, 4, 5]) expect(symmetric(n).derivedSubgroup().equals(alternating(n))).toBe(true);
    for (const p of [3, 5]) expect(GL(2, p).derivedSubgroup().equals(SL(2, p))).toBe(true);
    expect(alternating(5).isPerfect).toBe(true);
  });
  test("solvable, nilpotent and series", () => {
    expect(symmetric(4).isSolvable).toBe(true);
    expect(symmetric(5).isSolvable).toBe(false);
    expect(symmetric(3).isNilpotent).toBe(false);
    expect(heisenberg(3).nilpotencyClass).toBe(2);
    expect(dihedral(8).nilpotencyClass).toBe(3);
    expect(dihedral(6).nilpotencyClass).toBe(null);
    expect(symmetric(4).derivedSeries().map((H) => H.order)).toEqual([24, 12, 4, 1]);
    expect(dihedral(8).upperCentralSeries().map((H) => H.order)).toEqual([1, 2, 4, 16]);
    // every abelian group is nilpotent, every nilpotent group solvable
    for (const G of SMALL()) {
      if (G.isAbelian) expect(G.isNilpotent).toBe(true);
      if (G.isNilpotent) expect(G.isSolvable).toBe(true);
      if (G.isPGroup) expect(G.isNilpotent).toBe(true);
    }
  });
  test("composition series", () => {
    const S4 = symmetric(4).compositionSeries();
    expect(S4.factors).toEqual([2, 3, 2, 2]);
    expect(symmetric(5).compositionSeries().factors).toEqual([2, 60]);
    expect(alternating(5).isSimple).toBe(true);
    expect(PSL(2, 7).isSimple).toBe(true);
    expect(alternating(4).isSimple).toBe(false);
    expect(cyclic(7).isSimple).toBe(true);
    // Jordan–Hölder: the factors multiply to the order
    for (const G of SMALL()) expect(G.compositionSeries().factors.reduce((a, b) => a * b, 1)).toBe(G.order);
  });
});

describe("subgroups and quotients", () => {
  test("every listed subgroup is closed, distinct, and has order dividing |G| (Lagrange)", () => {
    for (const G of [...SMALL(), alternating(5)]) {
      const subs = G.subgroups();
      const ids = new Set();
      for (const H of subs) {
        expect(G.order % H.order).toBe(0);
        expect(isClosed(G, H.elements)).toBe(true);
        ids.add([...keys(G, H.elements)].sort().join("|"));
      }
      expect(ids.size).toBe(subs.length);
    }
  });
  test("known numbers of subgroups", () => {
    // Cₙ has one subgroup per divisor; Dₙ has τ(n) + σ(n)
    for (let n = 1; n <= 12; n += 1) {
      expect(cyclic(n).subgroups().length).toBe(divisors(n).length);
      expect(dihedral(n).subgroups().length).toBe(divisors(n).length + sigma(n));
    }
    expect(quaternion().subgroups().length).toBe(6);
    expect(symmetric(4).subgroups().length).toBe(30);
    expect(alternating(4).subgroups().length).toBe(10);
    expect(alternating(5).subgroups().length).toBe(59);
    expect(symmetric(5).subgroups().length).toBe(156);
    expect(directProduct(cyclic(2), cyclic(2), cyclic(2)).subgroups().length).toBe(16);
    expect(directProduct(cyclic(2), cyclic(2), cyclic(2), cyclic(2)).subgroups().length).toBe(67);
    expect(symmetric(4).subgroupClasses().length).toBe(11);
    expect(alternating(5).subgroupClasses().length).toBe(9);
    expect(PSL(2, 7).subgroupClasses().length).toBe(15);
  });
  test("normal subgroups against brute force", () => {
    for (const G of [...SMALL(), alternating(5), symmetric(5)]) {
      const brute = G.subgroups().filter((H) => bruteIsNormal(G, H));
      const fast = G.normalSubgroups();
      expect(fast.length).toBe(brute.length);
      for (const N of fast) {
        expect(bruteIsNormal(G, N)).toBe(true);
        expect(G.isNormal(N)).toBe(true);
      }
      for (const H of G.subgroups()) expect(G.isNormal(H)).toBe(bruteIsNormal(G, H));
    }
    expect(symmetric(4).normalSubgroups().map((N) => N.order)).toEqual([1, 4, 12, 24]);
    expect(quaternion().normalSubgroups().length).toBe(6); // every subgroup of Q₈ is normal
  });
  test("normalisers against brute force", () => {
    const G = symmetric(4);
    for (const H of G.subgroups()) {
      const hk = keys(G, H.elements);
      const brute = G.elements.filter((g) => H.elements.every((h) => hk.has(G.key(G.mul(G.mul(g, h), G.inv(g))))));
      expect(sameSet(G, G.normaliser(H).elements, brute)).toBe(true);
    }
  });
  test("cosets partition the group and quotients multiply cosets correctly", () => {
    for (const G of [dihedral(4), symmetric(4), quaternion(), alternating(4), heisenberg(3)]) {
      for (const H of G.subgroups()) {
        const cs = G.cosets(H);
        expect(cs.length * H.order).toBe(G.order);
        expect(new Set(cs.flatMap((c) => c.elements.map((x) => G.key(x)))).size).toBe(G.order);
      }
      for (const N of G.normalSubgroups()) {
        const Q = G.quotient(N);
        expect(Q.order).toBe(G.order / N.order);
        for (const a of G.elements) for (const b of G.elements) expect(Q.key(Q.mul(Q.cosetOf(a), Q.cosetOf(b)))).toBe(Q.key(Q.cosetOf(G.mul(a, b))));
        // left and right cosets agree for a normal subgroup
        const L = G.cosets(N).map((c) => [...keys(G, c.elements)].sort().join("|")).sort();
        const R = G.cosets(N, { side: "right" }).map((c) => [...keys(G, c.elements)].sort().join("|")).sort();
        expect(L).toEqual(R);
      }
    }
    const S3 = symmetric(3);
    const H = S3.subgroup(["(1 2)"]);
    expect(() => S3.quotient(H)).toThrow(/not a normal subgroup/);
  });
  test("intersections, joins and products", () => {
    const G = symmetric(4);
    const H = G.subgroup(["(1 2 3 4)"]);
    const K = G.subgroup(["(1 2)(3 4)", "(1 3)(2 4)"]);
    expect(G.intersection(H, K).order).toBe(2);
    expect(G.join(H, K).order).toBe(8);
    expect(G.product(H, K).length).toBe((H.order * K.order) / 2); // |HK| = |H||K|/|H ∩ K|
  });
});

describe("Sylow subgroups", () => {
  const primes = (n) => [2, 3, 5, 7, 11, 13].filter((p) => n % p === 0);
  test("Sylow's theorems on many groups, against counting subgroups", () => {
    for (const G of [...SMALL(), ...MEDIUM(), PSL(2, 7), dihedral(12), directProduct(symmetric(3), cyclic(3))]) {
      for (const p of primes(G.order)) {
        let pa = 1;
        while (G.order % (pa * p) === 0) pa *= p;
        const P = G.sylowSubgroup(p);
        expect(P.order).toBe(pa);
        expect(isClosed(G, P.elements)).toBe(true);
        const np = G.sylowCount(p);
        expect(np % p).toBe(1);
        expect((G.order / pa) % np).toBe(0);
        const all = G.sylowSubgroups(p);
        expect(all.length).toBe(np);
        if (G.order <= 168) expect(G.subgroups().filter((H) => H.order === pa).length).toBe(np);
      }
    }
  });
  test("known Sylow counts", () => {
    expect([2, 3, 5].map((p) => alternating(5).sylowCount(p))).toEqual([5, 10, 6]);
    expect([2, 3].map((p) => symmetric(4).sylowCount(p))).toEqual([3, 4]);
    expect([2, 3, 7].map((p) => PSL(2, 7).sylowCount(p))).toEqual([21, 28, 8]);
    // a group of order 15 has normal Sylow subgroups, so it is cyclic
    expect(cyclic(15).sylowCount(3)).toBe(1);
    // GL(2,3): its Sylow 2-subgroups are semidihedral of order 16
    expect(GL(2, 3).sylowCount(2)).toBe(3);
  });
});

describe("against algebraic's PermGroup, an independent implementation", () => {
  test("orders, subgroup counts, classes of subgroups, derived series and solvability agree", () => {
    const cases = [PermGroup.symmetric(4), PermGroup.alternating(5), PermGroup.from(["(1 2 3 4 5 6)", "(1 6)(2 5)(3 4)"], 6), PermGroup.from(["(1 2 3 4 5)", "(2 3 5 4)"], 5)];
    for (const P of cases) {
      const G = fromPermGroup(P);
      expect(G.order).toBe(P.order);
      expect(G.subgroups().length).toBe(P.subgroups().length);
      expect(G.subgroupClasses().map((c) => [c.order, c.conjugates]).sort()).toEqual(P.subgroupClasses().map((c) => [c.order, c.conjugates]).sort());
      expect(G.derivedSeries().map((H) => H.order)).toEqual(P.derivedSeries().map((H) => H.order));
      expect(G.isSolvable).toBe(P.isSolvable);
    }
  });
});

describe("products", () => {
  test("direct and semidirect products", () => {
    const P = directProduct(symmetric(3), cyclic(4));
    expect(P.order).toBe(24);
    expect(P.centre().order).toBe(4);
    // C₇ ⋊ C₃ with a ↦ a², the nonabelian group of order 21
    const F21 = metacyclic(7, 3, 2);
    expect(F21.order).toBe(21);
    expect(F21.isAbelian).toBe(false);
    expect(F21.sylowCount(3)).toBe(7);
    expect(() => metacyclic(7, 3, 3)).toThrow(/not an action/);
    // a map that is not an automorphism is refused
    expect(() => semidirectProduct(cyclic(4), cyclic(2), (h) => (h ? (x) => (x === 0 ? 0 : x === 1 ? 2 : x === 2 ? 1 : 3) : (x) => x))).toThrow(/homomorphism|one-to-one/);
    expect(matrixGroup([[[1, 1], [0, 1]]], 5).order).toBe(5);
  });
});
