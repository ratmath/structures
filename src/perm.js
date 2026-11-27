// The sign of a single permutation. Permutations are arrays of 0-based images
// (as fromCycles returns them); text such as "(1 2 3 4)" is read with fromCycles.

import { fromCycles } from "@ratmath/algebraic";

function readPerm(p, n) {
  if (typeof p === "string") {
    const size = n ?? Math.max(0, ...[...p.matchAll(/\d+/g)].map((m) => Number(m[0])));
    return fromCycles(p, size);
  }
  if (!Array.isArray(p)) throw new TypeError("a permutation is an array of images (0-based) or cycle text such as \"(1 2 3)\"");
  const seen = new Array(p.length).fill(false);
  for (const x of p) {
    if (!Number.isInteger(x) || x < 0 || x >= p.length || seen[x]) throw new RangeError(`[${p.join(", ")}] is not a permutation of 0..${p.length - 1}`);
    seen[x] = true;
  }
  return p;
}

/**
 * The sign of a permutation: 1 when it is a product of an even number of
 * transpositions, −1 when odd. A k-cycle is k − 1 transpositions, so the sign
 * is (−1)^(n − number of cycles).
 * @param {number[]|string} p  images 0..n−1, or cycle text with n
 * @param {number} [n]  points, for cycle text (default: the largest point named)
 * @returns {1|-1}
 */
export function sign(p, n) {
  const perm = readPerm(p, n);
  const seen = new Array(perm.length).fill(false);
  let cycles = 0;
  for (let i = 0; i < perm.length; i += 1) {
    if (seen[i]) continue;
    cycles += 1;
    for (let j = i; !seen[j]; j = perm[j]) seen[j] = true;
  }
  return (perm.length - cycles) % 2 === 0 ? 1 : -1;
}

/** True when the permutation is even (sign 1). */
export const isEven = (p, n) => sign(p, n) === 1;

/** True when the permutation is odd (sign −1). */
export const isOdd = (p, n) => sign(p, n) === -1;

/** The number of inversions: pairs i < j with p(i) > p(j). Its parity is the sign's. */
export function inversions(p, n) {
  const perm = readPerm(p, n);
  let count = 0;
  for (let i = 0; i < perm.length; i += 1) for (let j = i + 1; j < perm.length; j += 1) if (perm[i] > perm[j]) count += 1;
  return count;
}
