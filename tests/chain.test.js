import { describe, expect, test } from "bun:test";
import { chainComplex, ChainComplex, abelianGroup } from "../index.js";

const str = (C) => C.homologyGroups().map((H) => H.toString({ additive: true }));

// boundary matrices of a simplicial complex, built independently of the library
function simplicial(maxSimplices) {
  const faces = new Map();
  const key = (s) => s.join(",");
  for (const top of maxSimplices) {
    const s = [...top].sort((a, b) => a - b);
    for (let m = 1; m < 1 << s.length; m++) {
      const f = s.filter((_, i) => (m >> i) & 1);
      const d = f.length - 1;
      if (!faces.has(d)) faces.set(d, new Map());
      faces.get(d).set(key(f), f);
    }
  }
  const dims = [...faces.keys()].sort((a, b) => a - b);
  const lists = dims.map((d) => [...faces.get(d).values()]);
  const index = lists.map((l) => new Map(l.map((f, i) => [key(f), i])));
  const mats = [];
  for (let d = 1; d < lists.length; d++) {
    const M = lists[d - 1].map(() => Array(lists[d].length).fill(0));
    lists[d].forEach((f, j) => f.forEach((_, i) => {
      const g = f.filter((_, t) => t !== i);
      M[index[d - 1].get(key(g))][j] += i % 2 ? -1 : 1;
    }));
    mats.push(M);
  }
  return { mats, sizes: lists.map((l) => l.length) };
}

// rank over Q by exact fraction-free elimination (independent of Smith)
function rankQ(M) {
  const A = M.map((r) => r.map(BigInt));
  let r = 0;
  const cols = A.length ? A[0].length : 0;
  for (let c = 0; c < cols && r < A.length; c++) {
    const p = A.findIndex((row, i) => i >= r && row[c] !== 0n);
    if (p < 0) continue;
    [A[r], A[p]] = [A[p], A[r]];
    for (let i = r + 1; i < A.length; i++) {
      const f = A[i][c], g = A[r][c];
      A[i] = A[i].map((x, j) => x * g - A[r][j] * f);
    }
    r++;
  }
  return r;
}

describe("chainComplex: cell complexes from the textbook", () => {
  test("torus: Z, Z², Z", () => {
    const T = chainComplex([[[0, 0]], [[0], [0]]]);
    expect(str(T)).toEqual(["Z", "Z^2", "Z"]);
    expect(T.eulerCharacteristic).toBe(0);
  });
  test("projective plane: Z, Z/2, 0; over Z/2 every group is Z/2", () => {
    const P = chainComplex([[[0]], [[2]]]);
    expect(str(P)).toEqual(["Z", "Z/2", "0"]);
    expect(P.betti()).toEqual([1, 0, 0]);
    expect(P.betti({ coefficients: "Z/2" })).toEqual([1, 1, 1]);
    expect(P.betti({ coefficients: "Z/3" })).toEqual([1, 0, 0]);
    expect(String(P.homology(1, { coefficients: "Z/2" }))).toBe("Z/2");
  });
  test("Klein bottle: Z, Z ⊕ Z/2, 0", () => {
    const K = chainComplex([[[0, 0]], [[0], [2]]]);
    expect(str(K)).toEqual(["Z", "Z/2 ⊕ Z", "0"]);
    expect(K.betti({ coefficients: "Q" })).toEqual([1, 1, 0]);
  });
  test("torsion from ∂₂ with Smith form diag(2, 4)", () => {
    const C = chainComplex([[[0, 0]], [[2, 6], [4, 8]]]);
    expect(C.homology(1).toString()).toBe("Z/2 × Z/4");
    expect(C.invariantFactors(2)).toEqual([2n, 4n]);
  });
  test("real projective spaces RPⁿ: ∂ alternates 0 and 2", () => {
    for (let n = 1; n <= 6; n++) {
      const mats = Array.from({ length: n }, (_, i) => [[(i + 1) % 2 === 0 ? 2 : 0]]);
      const H = str(chainComplex(mats));
      for (let k = 1; k <= n; k++) {
        const want = k === n ? (n % 2 ? "Z" : "0") : k % 2 ? "Z/2" : "0";
        expect(H[k]).toBe(want);
      }
    }
  });
  test("refuses ∂∂ ≠ 0, mismatched sizes, and bad coefficients", () => {
    expect(() => chainComplex([[[1, 1]], [[1], [0]]])).toThrow(/not a chain complex/);
    expect(() => chainComplex([[[0, 0]], [[0], [0], [0]]])).toThrow();
    expect(() => chainComplex([[[0]]]).homology(0, { coefficients: "Z/4" })).toThrow(/prime/);
  });
  test("a point and empty maps with sizes", () => {
    expect(String(chainComplex([], { sizes: [1] }))).toBe("H0 = Z");
    const W = chainComplex([[], []], { sizes: [1, 0, 1] });  // S² as one 0-cell and one 2-cell
    expect(str(W)).toEqual(["Z", "0", "Z"]);
    expect(W.reducedHomology(0).toString({ additive: true })).toBe("0");
  });
});

describe("independent checks on simplicial complexes", () => {
  const named = {
    "S²": { tops: [[1, 2, 3], [1, 2, 4], [1, 3, 4], [2, 3, 4]], H: ["Z", "0", "Z"] },
    "S³": { tops: [[1, 2, 3, 4], [1, 2, 3, 5], [1, 2, 4, 5], [1, 3, 4, 5], [2, 3, 4, 5]], H: ["Z", "0", "0", "Z"] },
    "RP² (6 vertices)": { tops: [[1, 2, 3], [1, 3, 4], [1, 4, 5], [1, 5, 6], [1, 6, 2], [2, 3, 5], [3, 4, 6], [4, 5, 2], [5, 6, 3], [6, 2, 4]], H: ["Z", "Z/2", "0"] },
    "torus (7 vertices)": { tops: [[1, 2, 4], [2, 3, 5], [3, 4, 6], [4, 5, 7], [5, 6, 1], [6, 7, 2], [7, 1, 3], [1, 2, 6], [2, 3, 7], [3, 4, 1], [4, 5, 2], [5, 6, 3], [6, 7, 4], [7, 1, 5]], H: ["Z", "Z^2", "Z"] },
  };
  for (const [name, { tops, H }] of Object.entries(named)) {
    if (!H) continue;
    test(`${name}: known homology, χ and Betti numbers agree, universal coefficients hold`, () => {
      const { mats, sizes } = simplicial(tops);
      const C = chainComplex(mats);
      expect(C.sizes).toEqual(sizes);
      expect(str(C)).toEqual(H);
      // χ = alternating sum of Betti numbers over Q, Z/2 and Z/3
      for (const coefficients of ["Q", "Z/2", "Z/3"]) {
        const b = C.betti({ coefficients });
        expect(b.reduce((s, x, k) => s + (k % 2 ? -x : x), 0)).toBe(C.eulerCharacteristic);
      }
      // ranks over Q against independent elimination
      mats.forEach((M, i) => expect(C.rank(i + 1)).toBe(rankQ(M)));
      // dim H_k(Z/p) = b_k + t_k(p) + t_(k−1)(p)
      const groups = C.homologyGroups();
      for (const p of [2n, 3n]) {
        const t = groups.map((G) => G.invariants.filter((d) => d % p === 0n).length);
        const bp = C.betti({ coefficients: `Z/${p}` });
        groups.forEach((G, k) => expect(bp[k]).toBe(G.rank + t[k] + (k ? t[k - 1] : 0)));
      }
    });
  }
  test("ChainComplex is the class and toString lists the groups", () => {
    const C = chainComplex([[[0]], [[2]]]);
    expect(C instanceof ChainComplex).toBe(true);
    expect(String(C)).toBe("H0 = Z, H1 = Z/2, H2 = 0");
    expect(C.torsion.map(String)).toEqual(["1", "Z/2", "1"]);
    expect(abelianGroup([], 0).toString({ additive: true })).toBe("0");
  });
});
