// The standard groups, and ways of building new groups from old.

import { PermGroup, fromCycles, cycleString } from "@ratmath/algebraic";
import { Group, sup, isPrime, defaultKey } from "./group.js";
import { action } from "./action.js";

const compose = (a, b) => b.map((x) => a[x]); // a∘b: first b, then a
const range = (n) => [...Array(n).keys()];
const posInt = (n, what) => {
  if (!Number.isInteger(n) || n < 1) throw new RangeError(`${what} must be a positive whole number, not ${n}`);
  return n;
};

/** Cₙ: the integers mod n under addition, written multiplicatively as powers of a. */
export function cyclic(n) {
  posInt(n, "the order of a cyclic group");
  const G = Group.generate({ generators: n > 1 ? [1] : [], identity: 0, mul: (a, b) => (a + b) % n, format: (k) => (k === 0 ? "e" : k === 1 ? "a" : `a${sup(k)}`), name: `C${n}` });
  G.natural = action(G, range(n), (k, x) => (x + k) % n, { name: `C${n} turning ${n} beads` });
  return G;
}

/**
 * Dₙ: the symmetries of a regular n-gon, of order 2n, as rᵏsᶠ with r a turn
 * by 1/n and s a reflection (srs = r⁻¹). `natural` is its action on the
 * vertices 0, …, n−1.
 */
export function dihedral(n) {
  posInt(n, "the number of sides");
  const mul = ([a, f], [b, g]) => [(((a + (f ? -b : b)) % n) + n) % n, f ^ g];
  const fmt = ([k, f]) => (k === 0 && !f ? "e" : `${k === 0 ? "" : k === 1 ? "r" : `r${sup(k)}`}${f ? "s" : ""}`);
  const G = Group.generate({ generators: n > 1 ? [[1, 0], [0, 1]] : [[0, 1]], identity: [0, 0], mul, format: fmt, name: `D${n}` });
  G.natural = action(G, range(n), ([k, f], x) => (((f ? -x : x) + k) % n + n) % n, { name: `D${n} on the vertices of a ${n}-gon` });
  return G;
}

/**
 * The dicyclic group of order 4n: aᵏxᶠ with a²ⁿ = 1, x² = aⁿ, xax⁻¹ = a⁻¹.
 * n = 2 is the quaternion group Q₈, and n a power of 2 gives the generalised
 * quaternion groups.
 */
export function dicyclic(n) {
  posInt(n, "n");
  const m = 2 * n;
  const md = (k) => ((k % m) + m) % m;
  const mul = ([p, f], [q, g]) => (!f ? [md(p + q), g] : !g ? [md(p - q), 1] : [md(p - q + n), 0]);
  const fmt = ([k, f]) => (k === 0 && !f ? "e" : `${k === 0 ? "" : k === 1 ? "a" : `a${sup(k)}`}${f ? "x" : ""}`);
  return Group.generate({ generators: [[1, 0], [0, 1]], identity: [0, 0], mul, format: fmt, name: n === 2 ? "Q8" : `Dic${n}` });
}

/** The quaternion group Q₈ = {±1, ±i, ±j, ±k}. */
export function quaternion() {
  const G = dicyclic(2);
  const names = { "0,0": "1", "1,0": "i", "2,0": "−1", "3,0": "−i", "0,1": "j", "1,1": "k", "2,1": "−j", "3,1": "−k" };
  const H = Group.generate({ generators: G.generators, identity: [0, 0], mul: (a, b) => G.mul(a, b), format: (x) => names[x.join()], name: "Q8", labels: { "-1": [2, 0], "-i": [3, 0], "-j": [2, 1], "-k": [3, 1] } });
  return H;
}

/** The Klein four-group C₂ × C₂. */
export const kleinFour = () => directProduct(cyclic(2), cyclic(2), { name: "V4" });

// ---- permutation groups -----------------------------------------------------------------------------

/**
 * A permutation group from generators: arrays p with p[i] the image of i
 * (counting from 0), or cycle notation counted from 1 such as "(1 2 3)(4 5)".
 * Products follow algebraic's PermGroup: mul(a, b) = a∘b, first b then a.
 * `natural` is the action on the points 1, …, n.
 */
export function permutationGroup(generators, n, { name = null, maxOrder } = {}) {
  if (generators instanceof PermGroup) return permutationGroup(generators.generators, generators.degree, { name, maxOrder });
  if (n === undefined) n = Math.max(...generators.map((g) => (Array.isArray(g) ? g.length : Math.max(0, ...[...String(g).matchAll(/\d+/g)].map((m) => Number(m[0]))))));
  posInt(n, "the degree");
  const gens = generators.map((g) => (Array.isArray(g) ? [...g] : fromCycles(g, n)));
  for (const g of gens) if (g.length !== n || [...g].sort((a, b) => a - b).some((x, i) => x !== i)) throw new RangeError(`[${g}] is not a permutation of ${n} points`);
  const G = Group.generate({ generators: gens, identity: range(n), mul: compose, key: (p) => (Array.isArray(p) ? p.join(",") : typeof p === "string" ? (p.startsWith("(") ? fromCycles(p, n).join(",") : p) : String(p)), format: cycleString, name, maxOrder });
  G.degree = n;
  G.natural = action(G, range(n).map((i) => i + 1), (p, x) => p[x - 1] + 1, { name: `${G} on ${n} points` });
  G.toPermGroup = () => new PermGroup(G.generators, n);
  return G;
}

/** Sₙ, all permutations of n points. */
export function symmetric(n) {
  posInt(n, "n");
  const gens = n === 1 ? [] : n === 2 ? ["(1 2)"] : ["(1 2)", `(${range(n).map((i) => i + 1).join(" ")})`];
  return permutationGroup(gens, n, { name: `S${n}` });
}

/** Aₙ, the even permutations of n points. */
export function alternating(n) {
  posInt(n, "n");
  const gens = range(n).slice(2).map((i) => `(1 2 ${i + 1})`);
  return permutationGroup(gens, n, { name: `A${n}` });
}

/** algebraic's PermGroup as a Group (the same as permutationGroup(P)). */
export const fromPermGroup = (P, options) => permutationGroup(P.generators, P.degree, options);

// ---- matrix groups over F_p -------------------------------------------------------------------------

function primitiveRoot(p) {
  for (let g = 1; g < p; g += 1) {
    let x = g;
    let k = 1;
    while (x !== 1) { x = (x * g) % p; k += 1; }
    if (k === p - 1) return g;
  }
  return 1;
}

const matMul = (p) => (A, B) => A.map((row) => B[0].map((_, j) => row.reduce((s, a, k) => (s + a * B[k][j]) % p, 0)));
const ident = (n) => range(n).map((i) => range(n).map((j) => (i === j ? 1 : 0)));
const fmtMat = (A) => `[${A.map((r) => r.join(" ")).join("; ")}]`;

function checkField(n, p) {
  posInt(n, "the matrix size");
  if (!isPrime(p)) throw new RangeError(`matrix groups here are over F_p for a prime p, and ${p} is not prime`);
  if (p > 46340) throw new RangeError("p is too large for exact arithmetic with JavaScript numbers");
}

function transvections(n) {
  const out = [];
  for (let i = 0; i < n; i += 1) for (let j = 0; j < n; j += 1) {
    if (i === j) continue;
    const T = ident(n);
    T[i][j] = 1;
    out.push(T);
  }
  return out;
}

/**
 * The group generated by some invertible n×n matrices over F_p (entries
 * whole numbers 0 … p−1). `natural` acts on the column vectors of F_pⁿ.
 */
export function matrixGroup(generators, p, { name = null, maxOrder } = {}) {
  const n = generators[0]?.length ?? 1;
  checkField(n, p);
  const gens = generators.map((A) => A.map((r) => r.map((x) => ((Number(x) % p) + p) % p)));
  const G = Group.generate({ generators: gens, identity: ident(n), mul: matMul(p), key: (A) => (Array.isArray(A) ? A.map((r) => r.join(",")).join(";") : String(A)), format: fmtMat, name, maxOrder });
  const vectors = [];
  for (let v = 0; v < p ** n; v += 1) vectors.push(range(n).map((i) => Math.floor(v / p ** i) % p));
  G.natural = action(G, vectors, (A, v) => A.map((row) => row.reduce((s, a, k) => (s + a * v[k]) % p, 0)), { name: `${G} on F${p}^${n}` });
  G.field = { p, n };
  return G;
}

/** GL(n, p): all invertible n×n matrices over F_p. */
export function GL(n, p) {
  checkField(n, p);
  const D = ident(n);
  D[0][0] = primitiveRoot(p);
  return matrixGroup([...transvections(n), ...(p > 2 ? [D] : [])], p, { name: `GL(${n},${p})` });
}

/** SL(n, p): the matrices of determinant 1. */
export function SL(n, p) {
  checkField(n, p);
  return matrixGroup(n === 1 ? [[[1]]] : transvections(n), p, { name: `SL(${n},${p})` });
}

/** PSL(n, p) = SL(n, p) / its centre (the scalar matrices in it). Elements are cosets. */
export function PSL(n, p) {
  const S = SL(n, p);
  const Z = S.centre();
  Z.name = "Z";
  return S.quotient(Z, { name: `PSL(${n},${p})` });
}

/** The Heisenberg group: upper unitriangular 3×3 matrices over F_p, of order p³. */
export function heisenberg(p) {
  checkField(3, p);
  return matrixGroup([[[1, 1, 0], [0, 1, 0], [0, 0, 1]], [[1, 0, 0], [0, 1, 1], [0, 0, 1]]], p, { name: `Heis(${p})` });
}

// ---- products ---------------------------------------------------------------------------------------

/** G₁ × G₂ × …: tuples multiplied entry by entry. */
export function directProduct(...groups) {
  let opts = {};
  if (groups.length && !(groups[groups.length - 1] instanceof Group)) opts = groups.pop();
  if (!groups.length) throw new RangeError("a direct product needs at least one group");
  const k = groups.length;
  const gens = [];
  groups.forEach((G, i) => {
    for (const g of G.generators) gens.push(groups.map((H, j) => (j === i ? g : H.identity)));
  });
  const name = opts.name ?? (groups.every((G) => G.name) ? groups.map((G) => (/[×⋊]/.test(G.name) ? `(${G.name})` : G.name)).join(" × ") : null);
  const P = Group.generate({
    generators: gens,
    identity: groups.map((G) => G.identity),
    mul: (a, b) => groups.map((G, i) => G.mul(a[i], b[i])),
    key: (t) => (Array.isArray(t) && t.length === k ? t.map((x, i) => groups[i].key(x)).join("|") : defaultKey(t)),
    format: (t) => `(${t.map((x, i) => groups[i].format(x)).join(", ")})`,
    name,
  });
  P.factors = groups;
  return P;
}

/**
 * N ⋊ H, where `act(h)` gives the automorphism of N that h induces, as a
 * function on N's elements. Elements are pairs [n, h] with
 * (n₁, h₁)(n₂, h₂) = (n₁·act(h₁)(n₂), h₁h₂). Checks that act is a
 * homomorphism from H to the automorphisms of N.
 */
export function semidirectProduct(N, H, act, { name = null } = {}) {
  // tabulate each automorphism on N's element numbers
  const maps = H.elements.map((h) => {
    const f = act(h);
    const m = N.elements.map((x) => {
      const i = N.indexOf(f(x));
      if (i < 0) throw new RangeError(`act(${H.format(h)}) sends ${N.format(x)} outside N`);
      return i;
    });
    if (new Set(m).size !== N.order) throw new RangeError(`act(${H.format(h)}) is not one-to-one on N`);
    for (let a = 0; a < N.order; a += 1) for (const b of N.gens) {
      if (m[N.mulIndex(a, b)] !== N.mulIndex(m[a], m[b])) throw new RangeError(`act(${H.format(h)}) is not a homomorphism of N`);
    }
    return m;
  });
  for (let h1 = 0; h1 < H.order; h1 += 1) for (const h2 of H.gens) {
    const m12 = maps[H.mulIndex(h1, h2)];
    for (const x of N.gens) if (m12[x] !== maps[h1][maps[h2][x]]) throw new RangeError("act is not a homomorphism from H to Aut(N)");
  }
  const mul = ([n1, h1], [n2, h2]) => {
    const i1 = H.indexOf(h1);
    const n2i = maps[i1][N.indexOf(n2)];
    return [N.elements[N.mulIndex(N.indexOf(n1), n2i)], H.elements[H.mulIndex(i1, H.indexOf(h2))]];
  };
  const gens = [...N.generators.map((x) => [x, H.identity]), ...H.generators.map((h) => [N.identity, h])];
  return Group.generate({
    generators: gens,
    identity: [N.identity, H.identity],
    mul,
    key: (t) => (Array.isArray(t) && t.length === 2 ? `${N.key(t[0])}|${H.key(t[1])}` : defaultKey(t)),
    format: ([n, h]) => `(${N.format(n)}, ${H.format(h)})`,
    name: name ?? (N.name && H.name ? `${N.name} ⋊ ${H.name}` : null),
  });
}

/** Cₘ ⋊ Cₖ with the generator of Cₖ acting as x ↦ xʳ (needs rᵏ ≡ 1 mod m). */
export function metacyclic(m, k, r, { name = null } = {}) {
  let rk = 1;
  for (let i = 0; i < k; i += 1) rk = (rk * r) % m;
  if (rk !== 1 % m) throw new RangeError(`${r}^${k} is not 1 mod ${m}, so this is not an action`);
  const pw = (h) => { let t = 1; for (let i = 0; i < h; i += 1) t = (t * r) % m; return t; };
  const G = semidirectProduct(cyclic(m), cyclic(k), (h) => (x) => (x * pw(h)) % m, { name: name ?? `C${m} ⋊ C${k}` });
  return G;
}
