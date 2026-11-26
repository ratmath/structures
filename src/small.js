// The groups of order at most 20, one of each up to isomorphism, built from
// the named constructions. Their number for each order is the known count
// (OEIS A000001): 1, 1, 1, 2, 1, 2, 1, 5, 2, 2, 1, 5, 1, 2, 1, 14, 1, 5, 1, 5.
// The tests check that the ones listed for each order are pairwise
// non-isomorphic and that every group of those orders met elsewhere is one
// of them.

import { cyclic, dihedral, dicyclic, quaternion, alternating, symmetric, directProduct, semidirectProduct, metacyclic } from "./named.js";
import { isomorphism } from "./hom.js";

const C = cyclic;
const named = (G, name) => { G.name = name; return G; };
const x = (...gs) => directProduct(...gs);

function pauli() {
  // C4 ∘ D4: C4 × D4 with the central elements of order 2 identified
  const P = x(C(4), dihedral(4));
  const N = P.subgroup([[2, [2, 0]]]);
  N.name = "Z";
  return named(P.quotient(N), "C4 ∘ D4");
}

function c2c2SemiC4() {
  const V = x(C(2), C(2));
  return named(semidirectProduct(V, C(4), (h) => (h % 2 ? ([a, b]) => [b, a] : (v) => v)), "(C2 × C2) ⋊ C4");
}

function c3c3SemiC2() {
  const N = x(C(3), C(3));
  return named(semidirectProduct(N, C(2), (h) => (h ? ([a, b]) => [(3 - a) % 3, (3 - b) % 3] : (v) => v)), "(C3 × C3) ⋊ C2");
}

const TABLE = {
  1: [() => C(1)],
  2: [() => C(2)],
  3: [() => C(3)],
  4: [() => C(4), () => named(x(C(2), C(2)), "C2 × C2")],
  5: [() => C(5)],
  6: [() => C(6), () => named(symmetric(3), "S3")],
  7: [() => C(7)],
  8: [() => C(8), () => named(x(C(4), C(2)), "C4 × C2"), () => named(x(C(2), C(2), C(2)), "C2 × C2 × C2"), () => dihedral(4), () => quaternion()],
  9: [() => C(9), () => named(x(C(3), C(3)), "C3 × C3")],
  10: [() => C(10), () => dihedral(5)],
  11: [() => C(11)],
  12: [() => C(12), () => named(x(C(6), C(2)), "C6 × C2"), () => alternating(4), () => dihedral(6), () => dicyclic(3)],
  13: [() => C(13)],
  14: [() => C(14), () => dihedral(7)],
  15: [() => C(15)],
  16: [
    () => C(16),
    () => named(x(C(4), C(4)), "C4 × C4"),
    c2c2SemiC4,
    () => metacyclic(4, 4, 3, { name: "C4 ⋊ C4" }),
    () => named(x(C(8), C(2)), "C8 × C2"),
    () => metacyclic(8, 2, 5, { name: "M16" }),
    () => dihedral(8),
    () => metacyclic(8, 2, 3, { name: "SD16" }),
    () => named(dicyclic(4), "Q16"),
    () => named(x(C(4), C(2), C(2)), "C4 × C2 × C2"),
    () => named(x(C(2), dihedral(4)), "C2 × D4"),
    () => named(x(C(2), quaternion()), "C2 × Q8"),
    pauli,
    () => named(x(C(2), C(2), C(2), C(2)), "C2 × C2 × C2 × C2"),
  ],
  17: [() => C(17)],
  18: [() => C(18), () => named(x(C(6), C(3)), "C6 × C3"), () => dihedral(9), () => named(x(symmetric(3), C(3)), "S3 × C3"), c3c3SemiC2],
  19: [() => C(19)],
  20: [() => C(20), () => named(x(C(10), C(2)), "C10 × C2"), () => dihedral(10), () => dicyclic(5), () => metacyclic(5, 4, 2, { name: "C5 ⋊ C4" })],
};

/** Largest order the table covers. */
export const SMALL_ORDER_LIMIT = 20;

const cache = new Map();

/** One group of each isomorphism type of order n (n ≤ 20), each with its name. */
export function smallGroups(n) {
  if (!Number.isInteger(n) || n < 1 || n > SMALL_ORDER_LIMIT) throw new RangeError(`the table of small groups covers orders 1 to ${SMALL_ORDER_LIMIT}`);
  if (!cache.has(n)) cache.set(n, TABLE[n].map((make) => make()));
  return cache.get(n);
}

/** The number of groups of order n up to isomorphism, for n ≤ 20. */
export const numberOfGroups = (n) => smallGroups(n).length;

/**
 * Which group in the table G is isomorphic to: { name, group, iso }, or null
 * when |G| > 20.
 */
export function identify(G) {
  if (G.order > SMALL_ORDER_LIMIT) return null;
  for (const S of smallGroups(G.order)) {
    const iso = isomorphism(G, S);
    if (iso) return Object.freeze({ name: S.name, group: S, iso });
  }
  throw new Error(`internal: a group of order ${G.order} matches nothing in the table of small groups`);
}
