// Homomorphisms, isomorphism tests, automorphisms, the table of small
// groups, and group actions with orbit counting.
import { describe, expect, test } from "bun:test";
import {
  cyclic, dihedral, dicyclic, quaternion, kleinFour, symmetric, alternating, permutationGroup, fromPermGroup,
  GL, SL, PSL, heisenberg, directProduct, semidirectProduct, metacyclic,
  homomorphism, tryHomomorphism, isomorphism, isIsomorphic, automorphisms, automorphismGroup, innerAutomorphism, projection, inclusion,
  firstIsomorphism, secondIsomorphism, thirdIsomorphism,
  action, regularAction, conjugationAction, cosetAction, burnside,
  smallGroups, numberOfGroups, identify, PermGroup,
} from "../index.js";

const gcd = (a, b) => (b ? gcd(b, a % b) : a);
const phi = (n) => [...Array(n).keys()].filter((i) => gcd(i, n) === 1).length;
const divisors = (n) => [...Array(n).keys()].map((i) => i + 1).filter((d) => n % d === 0);

/** f is a bijective homomorphism, checked on every pair. */
function bruteIso(f) {
  const G = f.source;
  const H = f.target;
  if (new Set(G.elements.map((x) => H.key(f.apply(x)))).size !== H.order || G.order !== H.order) return false;
  return G.elements.every((a) => G.elements.every((b) => H.key(f.apply(G.mul(a, b))) === H.key(H.mul(f.apply(a), f.apply(b)))));
}

const parity = (p) => {
  let s = 0;
  const seen = new Array(p.length).fill(false);
  for (let i = 0; i < p.length; i += 1) {
    if (seen[i]) continue;
    let len = 0;
    for (let j = i; !seen[j]; j = p[j]) { seen[j] = true; len += 1; }
    s += len - 1;
  }
  return s % 2;
};

describe("homomorphisms", () => {
  test("the sign map Sₙ → C₂ and the first isomorphism theorem", () => {
    for (const n of [3, 4, 5]) {
      const S = symmetric(n);
      const sign = homomorphism(S, cyclic(2), (p) => parity(p));
      expect(sign.kernel().equals(alternating(n))).toBe(true);
      expect(sign.isSurjective).toBe(true);
      const { quotient, image, iso } = firstIsomorphism(sign);
      expect(quotient.order).toBe(2);
      expect(image.order).toBe(2);
      expect(bruteIso(iso)).toBe(true);
    }
  });
  test("the determinant GL(2, 3) → F₃ˣ has kernel SL(2, 3)", () => {
    const G = GL(2, 3);
    const det = homomorphism(G, cyclic(2), ([[a, b], [c, d]]) => (((a * d - b * c) % 3) + 3) % 3 === 1 ? 0 : 1);
    expect(det.kernel().equals(SL(2, 3))).toBe(true);
  });
  test("images of generators: well-defined or refused", () => {
    // a ↦ a² defines C₆ → C₃; a ↦ a does not define C₄ → C₃
    const f = homomorphism(cyclic(6), cyclic(3), [1]);
    expect(f.kernel().order).toBe(2);
    expect(tryHomomorphism(cyclic(4), cyclic(3), [1])).toBe(null);
    expect(() => homomorphism(cyclic(4), cyclic(3), [1])).toThrow(/do not define/);
    // the number of homomorphisms Cₘ → Cₙ is gcd(m, n)
    for (let m = 1; m <= 8; m += 1) for (let n = 1; n <= 8; n += 1) {
      const count = cyclic(n).elements.filter((y) => tryHomomorphism(cyclic(m), cyclic(n), m > 1 ? [y] : [])).length;
      expect(m > 1 ? count : 1).toBe(gcd(m, n));
    }
    // generators given as a Map
    const D = dihedral(4);
    const g = homomorphism(D, cyclic(2), new Map([["r", 0], ["s", 1]]));
    expect(g.kernel().order).toBe(4);
  });
  test("kernels are normal, images are subgroups, and homomorphisms compose", () => {
    const S4 = symmetric(4);
    // S4 acting on the three pairings of {1,2,3,4} gives S4 → S3 with kernel V4
    const pairings = [[[1, 2], [3, 4]], [[1, 3], [2, 4]], [[1, 4], [2, 3]]];
    const keyOf = (P) => P.map((p) => [...p].sort().join("")).sort().join("|");
    const A = action(S4, pairings, (g, P) => P.map((pair) => pair.map((x) => g[x - 1] + 1)), { key: keyOf });
    const S3 = fromPermGroup(A.toPermGroup());
    const f = homomorphism(S4, S3, (g) => A.permutation(g));
    expect(f.kernel().order).toBe(4);
    expect(S4.isNormal(f.kernel())).toBe(true);
    expect(f.isSurjective).toBe(true);
    const sign3 = homomorphism(S3, cyclic(2), (p) => parity(p));
    const h = f.then(sign3);
    expect(h.kernel().equals(alternating(4))).toBe(true);
    expect(f.preimage(S3.subgroup(["(1 2 3)"])).order).toBe(12);
  });
  test("projections and inclusions", () => {
    const G = dihedral(6);
    const N = G.subgroup(["r^2"]);
    const Q = G.quotient(N);
    const p = projection(Q);
    expect(p.kernel().equals(N)).toBe(true);
    expect(inclusion(N, G).isInjective).toBe(true);
  });
  test("second and third isomorphism theorems on examples", () => {
    const G = symmetric(4);
    const V = G.subgroup(["(1 2)(3 4)", "(1 3)(2 4)"]);
    const H = G.subgroup(["(1 2 3)", "(1 2)"]); // a copy of S3
    const two = secondIsomorphism(G, H, V);
    expect(two.HN.order).toBe(24);
    expect(two.HcapN.order).toBe(1);
    expect(bruteIso(two.iso)).toBe(true);
    const three = thirdIsomorphism(G, V, alternating(4));
    expect(three.left.order).toBe(2);
    expect(bruteIso(three.iso)).toBe(true);
  });
});

describe("isomorphism", () => {
  test("known isomorphisms, each checked on every pair of elements", () => {
    const pairs = [
      [symmetric(3), dihedral(3)], [GL(2, 2), symmetric(3)], [kleinFour(), dihedral(2)], [cyclic(6), directProduct(cyclic(2), cyclic(3))],
      [PSL(2, 5), alternating(5)], [PSL(2, 7), GL(3, 2)], [PSL(2, 3), alternating(4)], [quaternion(), dicyclic(2)],
      [SL(2, 2), symmetric(3)], [directProduct(cyclic(3), cyclic(5)), cyclic(15)], [metacyclic(8, 2, 7), dihedral(8)],
    ];
    for (const [G, H] of pairs) {
      const f = isomorphism(G, H);
      expect(f).not.toBe(null);
      expect(bruteIso(f)).toBe(true);
    }
  });
  test("known non-isomorphisms", () => {
    const pairs = [
      [cyclic(4), kleinFour()], [dihedral(4), quaternion()], [alternating(4), dihedral(6)], [SL(2, 3), symmetric(4)],
      [directProduct(cyclic(2), cyclic(4)), cyclic(8)], [dicyclic(3), dihedral(6)], [directProduct(cyclic(2), quaternion()), directProduct(cyclic(2), dihedral(4))],
      [heisenberg(3), metacyclic(9, 3, 4)], [SL(2, 5), symmetric(5)],
    ];
    for (const [G, H] of pairs) expect(isIsomorphic(G, H)).toBe(false);
  });
});

describe("the groups of order at most 20", () => {
  const A000001 = [1, 1, 1, 2, 1, 2, 1, 5, 2, 2, 1, 5, 1, 2, 1, 14, 1, 5, 1, 5];
  test("the right number of each order, pairwise non-isomorphic", () => {
    for (let n = 1; n <= 20; n += 1) {
      const gs = smallGroups(n);
      expect(gs.length).toBe(A000001[n - 1]);
      expect(numberOfGroups(n)).toBe(A000001[n - 1]);
      for (const G of gs) expect(G.order).toBe(n);
      for (let i = 0; i < gs.length; i += 1) for (let j = i + 1; j < gs.length; j += 1) expect(isIsomorphic(gs[i], gs[j])).toBe(false);
    }
  });
  test("every group of order ≤ 20 met inside larger groups is in the table", () => {
    // subgroups of S5, GL(2,3), PSL(2,7) and D8 × C2 found independently
    let checked = 0;
    for (const big of [symmetric(5), GL(2, 3), PSL(2, 7), directProduct(dihedral(8), cyclic(2)), directProduct(quaternion(), cyclic(2))]) {
      for (const c of big.subgroupClasses()) {
        if (c.order > 20) continue;
        const id = identify(c.group);
        expect(bruteIso(id.iso)).toBe(true);
        checked += 1;
      }
    }
    expect(checked).toBeGreaterThan(40);
    expect(identify(GL(2, 3).sylowSubgroup(2)).name).toBe("SD16");
    expect(identify(directProduct(dihedral(8), cyclic(2)).sylowSubgroup(2)) === null).toBe(true); // order 32 is beyond the table
    expect(identify(PSL(2, 3)).name).toBe("A4");
  });
});

describe("automorphisms", () => {
  test("|Aut(Cₙ)| = φ(n)", () => {
    for (let n = 1; n <= 20; n += 1) expect(automorphisms(cyclic(n)).length).toBe(phi(n));
  });
  test("known automorphism groups", () => {
    const sizes = [[kleinFour(), 6], [quaternion(), 24], [dihedral(4), 8], [symmetric(3), 6], [directProduct(cyclic(2), cyclic(2), cyclic(2)), 168], [symmetric(4), 24], [alternating(5), 120], [dihedral(5), 20]];
    for (const [G, n] of sizes) expect(automorphisms(G).length).toBe(n);
    const A = automorphismGroup(directProduct(cyclic(2), cyclic(2), cyclic(2)));
    expect(isIsomorphic(A, GL(3, 2))).toBe(true);
    expect(isIsomorphic(automorphismGroup(quaternion()), symmetric(4))).toBe(true);
    // Inn(G) ≅ G/Z(G), normal in Aut(G)
    for (const G of [dihedral(4), symmetric(4), quaternion(), alternating(4)]) {
      const Aut = automorphismGroup(G);
      expect(Aut.inner.order).toBe(G.order / G.centre().order);
      expect(Aut.isNormal(Aut.inner)).toBe(true);
    }
    for (const f of automorphisms(dihedral(4))) expect(bruteIso(f)).toBe(true);
    expect(innerAutomorphism(symmetric(3), "(1 2)").apply([1, 2, 0])).toEqual([2, 0, 1]);
  });
});

describe("actions", () => {
  test("orbit–stabiliser and orbits partitioning the points", () => {
    const cases = [symmetric(4).natural, dihedral(6).natural, GL(2, 3).natural, conjugationAction(symmetric(4)), cosetAction(alternating(5), alternating(5).sylowSubgroup(5))];
    for (const A of cases) {
      const orbs = A.orbits();
      expect(orbs.reduce((s, o) => s + o.length, 0)).toBe(A.degree);
      for (const x of A.points) expect(A.orbit(x).length * A.stabiliser(x).order).toBe(A.group.order);
      const c = A.countOrbits();
      expect(c.count).toBe(orbs.length);
    }
    // GL(2,3) on F₃²: the zero vector, and everything else
    expect(GL(2, 3).natural.orbits().map((o) => o.length).sort()).toEqual([1, 8]);
    // conjugation orbits are the conjugacy classes
    expect(conjugationAction(symmetric(4)).orbits().length).toBe(5);
  });
  test("regular and coset actions", () => {
    const G = dihedral(5);
    const R = regularAction(G);
    expect(R.isTransitive && R.isFree && R.isFaithful).toBe(true);
    // the kernel of the action on cosets of H is the core: the largest normal subgroup inside H
    const S4 = symmetric(4);
    const D = S4.sylowSubgroup(2);
    const C = cosetAction(S4, D);
    expect(C.degree).toBe(3);
    expect(C.kernel().order).toBe(4);
  });
  test("a function that is not an action is refused", () => {
    expect(() => action(cyclic(4), [0, 1, 2], (g, x) => (x + 1) % 3)).toThrow(/not an action/);
    expect(() => action(cyclic(2), [0, 1], (g, x) => (g ? 0 : x))).toThrow(/one-to-one|not an action/);
  });
  test("Burnside: necklaces and bracelets", () => {
    // necklaces: (1/n) Σ_{d|n} φ(d) k^(n/d); bracelets add the reflections
    const necklaces = (n, k) => divisors(n).reduce((s, d) => s + phi(d) * k ** (n / d), 0) / n;
    const bracelets = (n, k) => (n % 2 ? (necklaces(n, k) + k ** ((n + 1) / 2)) / 2 : (necklaces(n, k) + ((k + 1) * k ** (n / 2)) / 2) / 2);
    for (let n = 1; n <= 8; n += 1) for (let k = 1; k <= 3; k += 1) {
      expect(cyclic(n).natural.onColourings(k).countOrbits().count).toBe(necklaces(n, k));
      expect(dihedral(n).natural.onColourings(k).countOrbits().count).toBe(bracelets(n, k));
    }
  });
  test("Burnside: colouring the faces of a cube", () => {
    // the rotations of a cube on its six faces (as in discrete's symmetry.cubeFaces)
    const cube = permutationGroup([[0, 2, 3, 4, 1, 5], [1, 5, 2, 0, 4, 3]], 6);
    expect(cube.order).toBe(24);
    expect(isIsomorphic(cube, symmetric(4))).toBe(true);
    const counts = [1, 2, 3].map((k) => burnside(cube.natural.onColourings(k)).count);
    expect(counts).toEqual([1, 10, 57]);
    // with two colours, how many ways to paint exactly three faces red: 2 (adjacent at a corner, or a band)
    expect(cube.natural.onSubsets(3).orbits().length).toBe(2);
    // the cycle index route (discrete's countColourings) needs this image as a PermGroup
    expect(cube.natural.toPermGroup()).toBeInstanceOf(PermGroup);
    expect(cube.natural.toPermGroup().order).toBe(24);
  });
  test("fixed points, kernels and the class-by-class sum", () => {
    const A = dihedral(4).natural;
    expect(A.fixedPoints("s")).toEqual([0, 2]);
    const c = A.countOrbits();
    expect(c.total).toBe(8n);
    expect(String(c)).toBe("(1/8)·(4 + 0 + 2·2 + 2·0 + 2·0) = 1");
  });
});

void [heisenberg, semidirectProduct];
