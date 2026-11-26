// Homomorphisms between finite groups, isomorphism tests, automorphism groups,
// and the isomorphism theorems on examples.

import { Group } from "./group.js";

/**
 * Builds the map from images of generators (as element numbers), or returns
 * null when they do not extend to a homomorphism. Walks G from the identity
 * by right multiplication by generators: f(x·s) = f(x)·f(s); a clash means
 * no homomorphism has these images. When every edge agrees the map is a
 * homomorphism, because each element is a word in the generators.
 * `gens` may be a subset of G's elements; then the map is on ⟨gens⟩ only.
 */
function extend(G, H, gens, images, { injective = false } = {}) {
  const f = new Int32Array(G.order).fill(-1);
  f[0] = 0;
  const used = injective ? new Uint8Array(H.order) : null;
  if (used) used[0] = 1;
  const list = [0];
  for (let i = 0; i < list.length; i += 1) {
    const x = list[i];
    for (let k = 0; k < gens.length; k += 1) {
      const y = G.mulIndex(x, gens[k]);
      const fy = H.mulIndex(f[x], images[k]);
      if (f[y] >= 0) { if (f[y] !== fy) return null; continue; }
      if (used) { if (used[fy]) return null; used[fy] = 1; }
      f[y] = fy;
      list.push(y);
    }
  }
  return { map: f, size: list.length };
}

export class Homomorphism {
  #map;
  constructor(source, target, map, { name = null } = {}) {
    this.source = source;
    this.target = target;
    this.#map = Int32Array.from(map);
    this.name = name;
  }

  /** f(x). */
  apply(x) {
    const i = this.source.indexOf(x);
    if (i < 0) throw new RangeError(`${x} is not in the source group`);
    return this.target.elements[this.#map[i]];
  }
  applyIndex(i) { return this.#map[i]; }
  /** The images of all elements, by element number. */
  get table() { return [...this.#map]; }

  /** ker f = { x : f(x) = e }, a normal subgroup of the source. */
  kernel() {
    const list = [];
    for (let i = 0; i < this.source.order; i += 1) if (this.#map[i] === 0) list.push(i);
    return this.source.subgroupFromIndices(list, null, "K");
  }

  /** im f, a subgroup of the target. */
  image() {
    return this.target.subgroupFromIndices([...new Set(this.#map)], null);
  }

  /** f⁻¹(K) for a subgroup K of the target. */
  preimage(K) {
    const m = this.target.maskOf(K);
    const list = [];
    for (let i = 0; i < this.source.order; i += 1) if (m[this.#map[i]]) list.push(i);
    return this.source.subgroupFromIndices(list);
  }

  get isInjective() { return new Set(this.#map).size === this.source.order; }
  get isSurjective() { return new Set(this.#map).size === this.target.order; }
  get isIsomorphism() { return this.source.order === this.target.order && this.isInjective; }

  /** g∘f: first this, then g. */
  then(g) {
    return new Homomorphism(this.source, g.target, Array.from(this.#map, (y) => {
      const i = g.source === this.target ? y : g.source.indexOf(this.target.elements[y]);
      if (i < 0) throw new RangeError("the maps do not compose: an image lies outside the second map's source");
      return g.applyIndex(i);
    }));
  }

  inverse() {
    if (!this.isIsomorphism) throw new RangeError("only an isomorphism has an inverse");
    const inv = new Int32Array(this.source.order);
    this.#map.forEach((y, x) => { inv[y] = x; });
    return new Homomorphism(this.target, this.source, inv);
  }

  /** Where the source's generators go, as text. */
  toString() {
    const G = this.source;
    const H = this.target;
    return `${this.name ? `${this.name}: ` : ""}${G} → ${H}, ${G.generators.map((g) => `${G.format(g)} ↦ ${H.format(this.apply(g))}`).join(", ")}`;
  }
}

/**
 * A homomorphism G → H, given by a function on elements, or by the images of
 * G's generators (a list in the order of G.generators, or a Map from
 * elements that generate G to their images). Throws if no homomorphism fits.
 */
export function homomorphism(G, H, f, options = {}) {
  const res = tryHomomorphism(G, H, f, options);
  if (!res) throw new RangeError("these images do not define a homomorphism (a relation of the source group fails in the target)");
  return res;
}

/** As homomorphism, but returns null instead of throwing. */
export function tryHomomorphism(G, H, f, { name = null } = {}) {
  let gens;
  let images;
  if (typeof f === "function") {
    const map = G.elements.map((x) => {
      const i = H.indexOf(f(x));
      if (i < 0) throw new RangeError(`f(${G.format(x)}) is not in ${H}`);
      return i;
    });
    // f is a homomorphism when f(x·s) = f(x)f(s) for every x and generator s
    for (let x = 0; x < G.order; x += 1) for (const s of G.gens) if (map[G.mulIndex(x, s)] !== H.mulIndex(map[x], map[s])) return null;
    return new Homomorphism(G, H, map, { name });
  }
  if (f instanceof Map) {
    gens = [...f.keys()].map((g) => G.indexOf(g));
    if (gens.some((i) => i < 0)) throw new RangeError("a key of the map is not in the source group");
    images = [...f.values()].map((h) => H.indexOf(h));
  } else {
    if (f.length !== G.gens.length) throw new RangeError(`give one image for each of the ${G.gens.length} generators ${G.generators.map((g) => G.format(g)).join(", ")}`);
    gens = [...G.gens];
    images = f.map((h) => H.indexOf(h));
  }
  if (images.some((i) => i < 0)) throw new RangeError("an image is not in the target group");
  const ext = extend(G, H, gens, images);
  if (!ext) return null;
  if (ext.size !== G.order) throw new RangeError("those elements do not generate the source group");
  return new Homomorphism(G, H, ext.map, { name });
}

/** The identity map and the trivial map. */
export const identityMap = (G) => new Homomorphism(G, G, [...Array(G.order).keys()], { name: "id" });
export const trivialMap = (G, H) => new Homomorphism(G, H, new Array(G.order).fill(0));
/** The inclusion of a subgroup. */
export const inclusion = (H, G) => new Homomorphism(H, G, H.elements.map((x) => G.indexOf(x)));
/** G → G/N, x ↦ xN. */
export const projection = (Q) => new Homomorphism(Q.over, Q, Q.over.elements.map((x) => Q.indexOf(Q.cosetOf(x))), { name: "π" });

// ---- isomorphism -------------------------------------------------------------------------------------

function profile(G) {
  // invariants an isomorphism preserves, per element: (order, centraliser order)
  const cls = G.conjugacyClasses();
  const sig = new Array(G.order);
  for (const c of cls) for (const i of c.indices) sig[i] = `${c.order}:${c.size}`;
  return sig;
}

function invariants(G) {
  const sig = profile(G);
  const counts = new Map();
  for (const s of sig) counts.set(s, (counts.get(s) ?? 0) + 1);
  return [G.order, G.isAbelian, G.centre().order, G.derivedSubgroup().order, [...counts].sort().join(";")].join("/");
}

/** A short generating list chosen so each generator has few possible images. */
function searchGenerators(G, sig) {
  const counts = new Map();
  for (const s of sig) counts.set(s, (counts.get(s) ?? 0) + 1);
  const gens = [];
  let inside = G.closeIndices([]);
  while (inside.list.length < G.order) {
    let best = null;
    for (let x = 1; x < G.order; x += 1) {
      if (inside.inside[x]) continue;
      const K = G.closeIndices([...gens, x]);
      const score = [K.list.length, -counts.get(sig[x])];
      if (!best || score[0] > best.score[0] || (score[0] === best.score[0] && score[1] > best.score[1])) best = { x, K, score };
    }
    gens.push(best.x);
    inside = best.K;
  }
  return gens;
}

function* isoSearch(G, H, limit = Infinity) {
  if (G.order !== H.order) return;
  const sg = profile(G);
  const sh = profile(H);
  const gens = searchGenerators(G, sg);
  const candidates = gens.map((g) => {
    const out = [];
    for (let y = 0; y < H.order; y += 1) if (sh[y] === sg[g]) out.push(y);
    return out;
  });
  const images = new Array(gens.length);
  let found = 0;
  function* rec(k) {
    if (k === gens.length) {
      const ext = extend(G, H, gens, images, { injective: true });
      if (ext && ext.size === G.order) { found += 1; yield ext.map; }
      return;
    }
    for (const y of candidates[k]) {
      images[k] = y;
      // prune: the partial map must already be an injective homomorphism on ⟨g₁ … g_k⟩
      if (k + 1 < gens.length && !extend(G, H, gens.slice(0, k + 1), images.slice(0, k + 1), { injective: true })) continue;
      yield* rec(k + 1);
      if (found >= limit) return;
    }
  }
  yield* rec(0);
}

/** An isomorphism G → H, or null when the groups are not isomorphic. */
export function isomorphism(G, H) {
  if (G.order !== H.order) return null;
  if (invariants(G) !== invariants(H)) return null;
  for (const map of isoSearch(G, H, 1)) {
    const f = new Homomorphism(G, H, map);
    return homomorphism(G, H, G.generators.map((g) => f.apply(g))); // re-checked from scratch
  }
  return null;
}

export const isIsomorphic = (G, H) => isomorphism(G, H) !== null;

/** Every automorphism of G, as Homomorphisms. Refuses more than `max` (default 50,000). */
export function automorphisms(G, { max = 50000 } = {}) {
  const out = [];
  for (const map of isoSearch(G, G)) {
    out.push(new Homomorphism(G, G, map));
    if (out.length > max) throw new RangeError(`more than ${max} automorphisms`);
  }
  return out;
}

/**
 * Aut(G) as a group under composition. Elements are automorphisms (as
 * Homomorphisms), multiplied as (αβ)(x) = α(β(x)). `inner` is the subgroup of
 * conjugations x ↦ gxg⁻¹, a normal subgroup isomorphic to G/Z(G).
 */
export function automorphismGroup(G, options) {
  const all = automorphisms(G, options);
  const key = (a) => [...G.gens].map((s) => a.applyIndex(s)).join(",");
  const byKey = new Map(all.map((a) => [key(a), a]));
  const mul = (a, b) => byKey.get([...G.gens].map((s) => a.applyIndex(b.applyIndex(s))).join(","));
  const fmt = (a) => (key(a) === [...G.gens].join(",") ? "id" : G.generators.map((g) => `${G.format(g)}↦${G.format(a.apply(g))}`).join(", "));
  const A = Group.fromElements(all, mul, { key: (a) => (a instanceof Homomorphism ? key(a) : String(a)), format: fmt, name: `Aut(${G})` });
  const innerKeys = new Set();
  for (let g = 0; g < G.order; g += 1) innerKeys.add([...G.gens].map((s) => G.mulIndex(G.mulIndex(g, s), G.invIndex(g))).join(","));
  A.inner = A.subgroupFromIndices([...innerKeys].map((k) => A.indexOf(byKey.get(k))), null, `Inn(${G})`);
  A.of = G;
  return A;
}

/** The automorphism x ↦ gxg⁻¹. */
export const innerAutomorphism = (G, g) => homomorphism(G, G, (x) => G.conjugate(x, g));

// ---- the isomorphism theorems, on examples ----------------------------------------------------------

/**
 * First isomorphism theorem: for f: G → H, the map G/ker f → im f,
 * x·ker f ↦ f(x), is well defined and an isomorphism. Returns
 * { kernel, quotient, image, iso }, with iso checked.
 */
export function firstIsomorphism(f) {
  const K = f.kernel();
  const Q = f.source.quotient(K);
  const I = f.image();
  const iso = homomorphism(Q, I, (c) => f.apply(f.source.elements[c.repIndex]));
  if (!iso.isIsomorphism) throw new Error("internal: the first isomorphism theorem failed");
  return Object.freeze({ kernel: K, quotient: Q, image: I, iso });
}

/**
 * Second isomorphism theorem: for a subgroup H and a normal subgroup N of G,
 * HN is a subgroup, H ∩ N is normal in H, and H/(H ∩ N) ≅ HN/N by
 * h(H ∩ N) ↦ hN. Returns { HN, HcapN, left, right, iso }.
 */
export function secondIsomorphism(G, H, N) {
  if (!G.isNormal(N)) throw new RangeError("N must be normal in G");
  const HN = G.join(H, N);
  if (HN.order !== G.product(H, N).length) throw new Error("internal: HN is not the product set");
  const I = G.intersection(H, N);
  const Hs = G.subgroupFromIndices(H.elements.map((x) => G.indexOf(x)));
  const left = Hs.quotient(Hs.subgroupFromIndices(I.elements.map((x) => Hs.indexOf(x)), null, "H∩N"));
  const right = HN.quotient(HN.subgroupFromIndices(N.elements.map((x) => HN.indexOf(x)), null, "N"));
  const iso = homomorphism(left, right, (c) => right.cosetOf(Hs.elements[c.repIndex]));
  if (!iso.isIsomorphism) throw new Error("internal: the second isomorphism theorem failed");
  return Object.freeze({ HN, HcapN: I, left, right, iso });
}

/**
 * Third isomorphism theorem: for normal subgroups N ⊆ K of G,
 * (G/N)/(K/N) ≅ G/K by (gN)(K/N) ↦ gK. Returns { GN, KN, left, right, iso }.
 */
export function thirdIsomorphism(G, N, K) {
  if (!G.isNormal(N) || !G.isNormal(K) || !K.contains(N)) throw new RangeError("need normal subgroups N ⊆ K of G");
  const GN = G.quotient(N);
  const KN = GN.subgroupFromIndices([...new Set(K.elements.map((x) => GN.indexOf(GN.cosetOf(x))))], null, "K/N");
  const left = GN.quotient(KN);
  const right = G.quotient(K);
  const iso = homomorphism(left, right, (c) => right.cosetOf(G.elements[GN.elements[c.repIndex].repIndex]));
  if (!iso.isIsomorphism) throw new Error("internal: the third isomorphism theorem failed");
  return Object.freeze({ GN, KN, left, right, iso });
}
