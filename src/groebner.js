// Gröbner bases by Buchberger's algorithm.
//
//   const R = ring("x, y, z", { order: "lex" });
//   const G = groebner([R.poly("x^2 + y^2 + z^2 - 1"), R.poly("x - y"), R.poly("y - z^2")]);
//   G.basis;             // the reduced Gröbner basis: unique for the ideal and the order
//   G.reduce(f);         // the normal form of f: zero exactly when f is in the ideal
//   G.express(f);        // f as a combination of the ORIGINAL generators
//   G.verify();          // replays the proof that this is a Gröbner basis of that ideal
//
// The pair selection is the normal strategy with the Gebauer–Möller
// criteria. Over Q the work is done on primitive integer polynomials, so no
// fractions are carried; with `certificates` (or the first time express is
// called) the cofactors of every basis element over the generators are
// carried along in field arithmetic, which is slower.

import { MPoly, PolyRing } from "./mpoly.js";
import { QQ } from "./field.js";
import { keyOf, dividesExp, divExp, mulExp, lcmExp, orderOf } from "./monomial.js";

const babs = (a) => (a < 0n ? -a : a);
const bgcd = (a, b) => {
  a = babs(a);
  b = babs(b);
  while (b) [a, b] = [b, a % b];
  return a;
};

/** Coefficient arithmetic: BigInt (fraction-free, for Q) or the ring's field. */
function domain(K, integer) {
  if (integer) return { integer: true, isZero: (a) => a === 0n, sub: (a, b) => a - b, mul: (a, b) => a * b, neg: (a) => -a, one: 1n };
  return { integer: false, isZero: (a) => K.sign(a) === 0, sub: K.sub, mul: K.mul, neg: K.neg, div: K.div, add: K.add, one: K.one };
}

const ZEROP = Object.freeze({ e: [], c: [] });

/** a·x^ma·f − b·x^mb·g for sorted term lists; a, b null mean 1 and ma, mb null mean no shift. */
function lin(D, cmp, a, ma, f, b, mb, g) {
  const E = [];
  const C = [];
  const nf = f.e.length;
  const ng = g.e.length;
  let i = 0;
  let j = 0;
  let ef = i < nf ? (ma ? mulExp(f.e[0], ma) : f.e[0]) : null;
  let eg = j < ng ? (mb ? mulExp(g.e[0], mb) : g.e[0]) : null;
  while (ef || eg) {
    const s = !eg ? 1 : !ef ? -1 : cmp(ef, eg);
    if (s > 0) {
      E.push(ef);
      C.push(a === null ? f.c[i] : D.mul(a, f.c[i]));
      i += 1;
      ef = i < nf ? (ma ? mulExp(f.e[i], ma) : f.e[i]) : null;
    } else if (s < 0) {
      E.push(eg);
      C.push(D.neg(b === null ? g.c[j] : D.mul(b, g.c[j])));
      j += 1;
      eg = j < ng ? (mb ? mulExp(g.e[j], mb) : g.e[j]) : null;
    } else {
      const v = D.sub(a === null ? f.c[i] : D.mul(a, f.c[i]), b === null ? g.c[j] : D.mul(b, g.c[j]));
      if (!D.isZero(v)) {
        E.push(ef);
        C.push(v);
      }
      i += 1;
      j += 1;
      ef = i < nf ? (ma ? mulExp(f.e[i], ma) : f.e[i]) : null;
      eg = j < ng ? (mb ? mulExp(g.e[j], mb) : g.e[j]) : null;
    }
  }
  return { e: E, c: C };
}

/** general product of two term lists (field arithmetic) */
function mulP(D, cmp, f, g) {
  if (!f.e.length || !g.e.length) return ZEROP;
  const m = new Map();
  for (let i = 0; i < f.e.length; i += 1)
    for (let j = 0; j < g.e.length; j += 1) {
      const e = mulExp(f.e[i], g.e[j]);
      const k = keyOf(e);
      const v = D.mul(f.c[i], g.c[j]);
      const t = m.get(k);
      if (t) t.c = D.add(t.c, v);
      else m.set(k, { e, c: v });
    }
  const ts = [...m.values()].filter((t) => !D.isZero(t.c)).sort((x, y) => cmp(y.e, x.e));
  return { e: ts.map((t) => t.e), c: ts.map((t) => t.c) };
}

function primitive(p) {
  if (!p.c.length) return p;
  let g = 0n;
  for (const c of p.c) {
    g = bgcd(g, c);
    if (g === 1n) break;
  }
  if (p.c[0] < 0n) g = -g;
  return g === 1n ? p : { e: p.e, c: p.c.map((c) => c / g) };
}
function monic(D, p) {
  if (!p.c.length) return { p, lead: null };
  const lead = p.c[0];
  return { p: { e: p.e, c: p.c.map((c) => D.div(c, lead)) }, lead };
}

/**
 * The normal form of f modulo the reducers (entries { p, lm }): no term of
 * the result is divisible by a leading monomial. In integer mode the result
 * is a multiple of the true remainder. With `quots`, the quotients are
 * collected (field mode): f = Σ quots[i]·reducers[i] + remainder.
 */
/** When set, the time (Date.now) after which a running computation gives up. */
let deadline = 0;
const checkTime = () => {
  if (deadline && Date.now() > deadline) throw new RangeError("the Gröbner basis computation passed its time limit");
};

function normalForm(D, cmp, f, reducers, quots = null) {
  let p = f;
  const RE = [];
  let RC = [];
  let steps = 0;
  while (p.e.length) {
    if (deadline) checkTime();
    const lm = p.e[0];
    let hit = -1;
    for (let k = 0; k < reducers.length; k += 1) {
      if (dividesExp(reducers[k].lm, lm)) {
        hit = k;
        break;
      }
    }
    if (hit < 0) {
      RE.push(lm);
      RC.push(p.c[0]);
      p = { e: p.e.slice(1), c: p.c.slice(1) };
      continue;
    }
    const g = reducers[hit].p;
    const m = divExp(lm, g.e[0]);
    if (D.integer) {
      const d = bgcd(p.c[0], g.c[0]);
      let a = g.c[0] / d;
      let b = p.c[0] / d;
      if (a < 0n) {
        a = -a;
        b = -b;
      }
      p = lin(D, cmp, a === 1n ? null : a, null, p, b, m, g);
      if (a !== 1n) RC = RC.map((c) => c * a);
      steps += 1;
      if (steps % 12 === 0) {
        // keep the numbers small: remove the common content of what is left and what is done
        let c0 = 0n;
        for (const c of p.c) {
          c0 = bgcd(c0, c);
          if (c0 === 1n) break;
        }
        if (c0 !== 1n) for (const c of RC) {
          c0 = bgcd(c0, c);
          if (c0 === 1n) break;
        }
        if (c0 > 1n) {
          p = { e: p.e, c: p.c.map((c) => c / c0) };
          RC = RC.map((c) => c / c0);
        }
      }
    } else {
      const q = D.div(p.c[0], g.c[0]);
      p = lin(D, cmp, null, null, p, q, m, g);
      if (quots) {
        const t = quots[hit] ?? (quots[hit] = { e: [], c: [] });
        t.e.push(m);
        t.c.push(q);
      }
    }
  }
  return { e: RE, c: RC };
}

const coprimeExp = (a, b) => {
  for (let i = 0; i < a.length; i += 1) if (a[i] && b[i]) return false;
  return true;
};
const sameExp = (a, b) => {
  for (let i = 0; i < a.length; i += 1) if (a[i] !== b[i]) return false;
  return true;
};

function toInternal(R, f, integer) {
  const ts = f._terms();
  if (!integer) return { e: ts.map((t) => t.e), c: ts.map((t) => t.c) };
  let L = 1n;
  for (const t of ts) L = (L * t.c[1]) / bgcd(L, t.c[1]);
  return primitive({ e: ts.map((t) => t.e), c: ts.map((t) => t.c[0] * (L / t.c[1])) });
}
function toPublic(R, p, integer) {
  const m = new Map();
  if (integer) {
    const lead = p.c[0];
    for (let i = 0; i < p.e.length; i += 1) {
      const g = bgcd(p.c[i], lead);
      let [n, d] = [p.c[i] / g, lead / g];
      if (d < 0n) [n, d] = [-n, -d];
      m.set(keyOf(p.e[i]), { e: p.e[i], c: [n, d] });
    }
  } else for (let i = 0; i < p.e.length; i += 1) m.set(keyOf(p.e[i]), { e: p.e[i], c: p.c[i] });
  return new MPoly(R, m);
}
const fieldPoly = (R, p) => {
  const m = new Map();
  for (let i = 0; i < p.e.length; i += 1) m.set(keyOf(p.e[i]), { e: p.e[i], c: p.c[i] });
  return new MPoly(R, m);
};

/**
 * Buchberger's algorithm. gens: internal polynomials. With `track`, each
 * element carries cof[k], its cofactor on generator k (field mode only).
 */
function buchberger(D, cmp, gens, options = {}) {
  const before = deadline;
  deadline = options.timeLimit ? Date.now() + options.timeLimit : 0;
  try {
    return buchbergerRun(D, cmp, gens, options);
  } finally {
    deadline = before;
  }
}

function buchbergerRun(D, cmp, gens, { track = false, maxPairs = 200000, trace = null } = {}) {
  const all = []; // { p, lm, alive, cof }
  let pairs = []; // { i, j, lcm }
  const stats = { pairs: 0, zero: 0, skipped: 0, added: 0 };
  const nGen = gens.length;
  const alive = () => all.filter((g) => g.alive);

  const combineCof = (base, quots, reducers) => {
    // base − Σ quots[i]·reducers[i].cof
    const out = base.slice();
    quots.forEach((q, i) => {
      if (!q) return;
      const rc = reducers[i].cof;
      for (let k = 0; k < nGen; k += 1) {
        if (!rc[k].e.length) continue;
        out[k] = lin(D, cmp, null, null, out[k], null, null, mulP(D, cmp, q, rc[k]));
      }
    });
    return out;
  };

  const add = (h, cof) => {
    const idx = all.length;
    const lm = h.e[0];
    // Gebauer–Möller update
    const C = [];
    all.forEach((g, i) => {
      if (g.alive) C.push({ i, j: idx, lcm: lcmExp(g.lm, lm), coprime: coprimeExp(g.lm, lm) });
    });
    const Dl = [];
    for (let a = 0; a < C.length; a += 1) {
      const c = C[a];
      let keep = c.coprime;
      if (!keep) {
        keep = true;
        for (let b = a + 1; b < C.length && keep; b += 1) if (dividesExp(C[b].lcm, c.lcm)) keep = false;
        for (let b = 0; b < Dl.length && keep; b += 1) if (dividesExp(Dl[b].lcm, c.lcm)) keep = false;
      }
      if (keep) Dl.push(c);
      else stats.skipped += 1;
    }
    const E = Dl.filter((c) => {
      if (c.coprime) stats.skipped += 1;
      return !c.coprime;
    });
    pairs = pairs.filter((pr) => {
      const ok = !dividesExp(lm, pr.lcm) || sameExp(lcmExp(all[pr.i].lm, lm), pr.lcm) || sameExp(lcmExp(all[pr.j].lm, lm), pr.lcm);
      if (!ok) stats.skipped += 1;
      return ok;
    });
    for (const c of E) pairs.push({ i: c.i, j: c.j, lcm: c.lcm });
    for (const g of all) if (g.alive && dividesExp(lm, g.lm)) g.alive = false;
    all.push({ p: h, lm, alive: true, cof });
    stats.added += 1;
    return idx;
  };

  const finish = (r, cofBase, quots, reducers) => {
    if (!r.e.length) return null;
    if (D.integer) return { p: primitive(r), cof: null };
    const { p, lead } = monic(D, r);
    let cof = null;
    if (track) {
      cof = combineCof(cofBase, quots, reducers).map((c) => (c.e.length ? { e: c.e, c: c.c.map((x) => D.div(x, lead)) } : c));
    }
    return { p, cof };
  };

  const unit = (k) => Array.from({ length: nGen }, (_, i) => (i === k ? { e: [new Array(gens[0].e[0].length).fill(0)], c: [D.one] } : ZEROP));

  // the generators, smallest leading monomial first
  const order = gens.map((g, k) => k).sort((a, b) => cmp(gens[a].e[0], gens[b].e[0]));
  for (const k of order) {
    const reducers = alive();
    const quots = track ? [] : null;
    const r = normalForm(D, cmp, gens[k], reducers, quots);
    const h = finish(r, track ? unit(k) : null, quots, reducers);
    const idx = h ? add(h.p, h.cof) : null;
    if (trace) trace.push({ generator: k, result: idx });
  }
  while (pairs.length) {
    // the normal strategy: the pair with the smallest lcm
    let best = 0;
    for (let k = 1; k < pairs.length; k += 1) if (cmp(pairs[k].lcm, pairs[best].lcm) < 0) best = k;
    const pr = pairs[best];
    pairs[best] = pairs[pairs.length - 1];
    pairs.pop();
    stats.pairs += 1;
    if (stats.pairs > maxPairs) throw new RangeError(`the Gröbner basis computation passed ${maxPairs} S-polynomials; the problem is too large for this order (try "grevlex")`);
    const [f, g] = [all[pr.i], all[pr.j]];
    const [mf, mg] = [divExp(pr.lcm, f.lm), divExp(pr.lcm, g.lm)];
    let s;
    let cofS = null;
    if (D.integer) {
      const d = bgcd(f.p.c[0], g.p.c[0]);
      s = lin(D, cmp, g.p.c[0] / d, mf, f.p, f.p.c[0] / d, mg, g.p);
    } else {
      s = lin(D, cmp, null, mf, f.p, null, mg, g.p); // both monic
      if (track) cofS = f.cof.map((c, k) => lin(D, cmp, null, mf, c, null, mg, g.cof[k]));
    }
    const reducers = alive();
    const quots = track ? [] : null;
    const r = normalForm(D, cmp, s, reducers, quots);
    const h = finish(r, cofS, quots, reducers);
    if (!h) stats.zero += 1;
    const idx = h ? add(h.p, h.cof) : null;
    if (trace) trace.push({ pair: [pr.i, pr.j], lcm: pr.lcm, result: idx });
  }

  // the reduced basis: minimal leading monomials, tails reduced
  let G = alive();
  G = G.filter((g, a) => !G.some((h, b) => b !== a && dividesExp(h.lm, g.lm) && (!sameExp(h.lm, g.lm) || b < a)));
  G.sort((a, b) => cmp(a.lm, b.lm));
  const out = G.map((g, a) => {
    const others = G.filter((_, b) => b !== a);
    const quots = track ? [] : null;
    const r = normalForm(D, cmp, g.p, others, quots);
    if (D.integer) return { p: primitive(r), lm: g.lm, cof: null };
    const { p, lead } = monic(D, r);
    const cof = track ? combineCof(g.cof, quots, others).map((c) => (c.e.length ? { e: c.e, c: c.c.map((x) => D.div(x, lead)) } : c)) : null;
    return { p, lm: g.lm, cof };
  });
  return { basis: out, history: all, stats };
}

/** The S-polynomial of two polynomials: the combination in which their leading terms cancel. */
export function sPolynomial(f, g) {
  const R = f.ring;
  g = R.from(g);
  if (f.isZero() || g.isZero()) throw new RangeError("the S-polynomial needs two non-zero polynomials");
  const [lf, lg] = [f._lead(), g._lead()];
  const L = lcmExp(lf.e, lg.e);
  const mono = (e, c) => new MPoly(R, new Map([[keyOf(e), { e, c }]]));
  return mono(divExp(L, lf.e), R.K.inv(lf.c)).mul(f).sub(mono(divExp(L, lg.e), R.K.inv(lg.c)).mul(g));
}

/**
 * A comparison function that eliminates the given variables: any monomial
 * containing one of them is larger than any monomial without. Inside each
 * block the order is graded reverse lexicographic.
 */
export function eliminationOrder(n, eliminate) {
  const first = eliminate;
  const rest = Array.from({ length: n }, (_, i) => i).filter((i) => !first.includes(i));
  const block = (idx) => (a, b) => {
    let d = 0;
    for (const i of idx) d += a[i] - b[i];
    if (d !== 0) return d;
    for (let k = idx.length - 1; k >= 0; k -= 1) {
      const i = idx[k];
      if (a[i] !== b[i]) return b[i] - a[i];
    }
    return 0;
  };
  const [c1, c2] = [block(first), block(rest)];
  return (a, b) => c1(a, b) || c2(a, b);
}

export class GroebnerBasis {
  /** Use groebner(polys, options). */
  constructor(R, generators, options) {
    this.ring = R;
    this.generators = generators;
    this._options = options;
    this._integer = R.K === QQ;
    const cmp = R.cmp;
    this._cmp = cmp;
    this._D = domain(R.K, false);
    const trace = [];
    if (!generators.length) {
      this._basis = [];
      this._certified = true;
      this.stats = { pairs: 0, zero: 0, skipped: 0, added: 0 };
      this.trace = trace;
    } else if (options.certificates || !this._integer) {
      const res = buchberger(this._D, cmp, generators.map((g) => toInternal(R, g, false)), { track: !!options.certificates, maxPairs: options.maxPairs, timeLimit: options.timeLimit, trace });
      this._basis = res.basis;
      this._certified = !!options.certificates;
      this.stats = res.stats;
      this.trace = trace;
    } else {
      const DZ = domain(R.K, true);
      const res = buchberger(DZ, cmp, generators.map((g) => toInternal(R, g, true)), { maxPairs: options.maxPairs, timeLimit: options.timeLimit, trace });
      // back to monic polynomials over Q
      this._basis = res.basis.map((g) => {
        const lead = g.p.c[0];
        return { lm: g.lm, cof: null, p: { e: g.p.e, c: g.p.c.map((c) => {
          const d = bgcd(c, lead);
          return lead < 0n ? [-c / d, -lead / d] : [c / d, lead / d];
        }) } };
      });
      this._certified = false;
      this.stats = res.stats;
      this.trace = trace;
    }
    /** The reduced Gröbner basis: monic, smallest leading monomial first. */
    this.basis = Object.freeze(this._basis.map((g) => fieldPoly(R, g.p)));
  }
  get length() {
    return this.basis.length;
  }
  /** The leading monomials of the basis, as exponent lists. */
  get leadingExponents() {
    return this._basis.map((g) => [...g.lm]);
  }
  get leadingMonomials() {
    return this.basis.map((g) => g.leadingMonomial());
  }
  /** True when the ideal is everything (the basis is {1}). */
  get isTrivial() {
    return this._basis.length === 1 && this._basis[0].lm.every((x) => x === 0);
  }
  get isZero() {
    return this._basis.length === 0;
  }
  _in(f) {
    return toInternal(this.ring, this.ring.from(f), false);
  }
  /** The normal form of f: the unique remainder on division by the basis. */
  reduce(f) {
    return fieldPoly(this.ring, normalForm(this._D, this._cmp, this._in(f), this._basis));
  }
  /** f = Σ quotients[i]·basis[i] + remainder. */
  divide(f) {
    const quots = [];
    const r = normalForm(this._D, this._cmp, this._in(f), this._basis, quots);
    return { quotients: this._basis.map((_, i) => (quots[i] ? fieldPoly(this.ring, quots[i]) : this.ring.zero)), remainder: fieldPoly(this.ring, r) };
  }
  /** Is f in the ideal? */
  contains(f) {
    return normalForm(this._D, this._cmp, this._in(f), this._basis).e.length === 0;
  }
  _certify() {
    if (this._certified) return;
    const res = buchberger(this._D, this._cmp, this.generators.map((g) => toInternal(this.ring, g, false)), { track: true, maxPairs: this._options.maxPairs, timeLimit: this._options.timeLimit });
    // the reduced basis is unique, so the certified run has the same elements in the same order
    this._basis = res.basis;
    this._certified = true;
  }
  /** Each basis element as a combination of the generators: cofactors[i][k] multiplies generators[k]. */
  get cofactors() {
    this._certify();
    return this._basis.map((g) => g.cof.map((c) => fieldPoly(this.ring, c)));
  }
  /**
   * f as a combination of the original generators:
   * { cofactors, check } with f = Σ cofactors[k]·generators[k]. Throws when
   * f is not in the ideal.
   */
  express(f) {
    const F = this.ring.from(f);
    this._certify();
    const quots = [];
    const r = normalForm(this._D, this._cmp, this._in(F), this._basis, quots);
    if (r.e.length) throw new RangeError("that polynomial is not in the ideal: its normal form is " + fieldPoly(this.ring, r).toString());
    const D = this._D;
    const cmp = this._cmp;
    const n = this.generators.length;
    let out = Array.from({ length: n }, () => ZEROP);
    quots.forEach((q, i) => {
      if (!q) return;
      const cof = this._basis[i].cof;
      out = out.map((o, k) => (cof[k].e.length ? lin(D, cmp, null, null, o, D.neg(D.one), null, mulP(D, cmp, q, cof[k])) : o));
    });
    const cofactors = out.map((c) => fieldPoly(this.ring, c));
    let sum = this.ring.zero;
    cofactors.forEach((c, k) => {
      sum = sum.add(c.mul(this.generators[k]));
    });
    return { cofactors, check: sum.equals(F) };
  }
  /**
   * Replays the proof. Three facts make this basis THE Gröbner basis of the
   * ideal of the generators: every generator reduces to zero (so the ideal
   * of the basis contains them), every S-polynomial of two basis elements
   * reduces to zero (Buchberger's criterion), and every basis element is a
   * combination of the generators (checked from the cofactors).
   */
  verify() {
    const generatorsReduce = this.generators.every((g) => this.contains(g));
    let sPairsReduce = true;
    let pairs = 0;
    for (let i = 0; i < this.basis.length && sPairsReduce; i += 1)
      for (let j = i + 1; j < this.basis.length && sPairsReduce; j += 1) {
        if (coprimeExp(this._basis[i].lm, this._basis[j].lm)) continue;
        pairs += 1;
        if (!this.contains(sPolynomial(this.basis[i], this.basis[j]))) sPairsReduce = false;
      }
    const cof = this.cofactors;
    const basisInIdeal = this.basis.every((g, i) => {
      let sum = this.ring.zero;
      cof[i].forEach((c, k) => {
        sum = sum.add(c.mul(this.generators[k]));
      });
      return sum.equals(g);
    });
    return { generatorsReduce, sPairsReduce, basisInIdeal, sPairsChecked: pairs, ok: generatorsReduce && sPairsReduce && basisInIdeal };
  }
  /** The largest number of variables that are independent modulo the ideal: its dimension (−1 for the whole ring). */
  get dimension() {
    if (this._dim !== undefined) return this._dim;
    if (this.isTrivial) return (this._dim = -1);
    const n = this.ring.n;
    const lms = this._basis.map((g) => g.lm);
    let best = 0;
    for (let mask = 0; mask < 1 << n; mask += 1) {
      let size = 0;
      for (let i = 0; i < n; i += 1) if (mask & (1 << i)) size += 1;
      if (size <= best) continue;
      // independent: no leading monomial uses only variables of the set
      if (lms.every((e) => e.some((x, i) => x > 0 && !(mask & (1 << i))))) best = size;
    }
    return (this._dim = best);
  }
  get isZeroDimensional() {
    return this.dimension === 0;
  }
  /**
   * The monomials that are not divisible by any leading monomial, smallest
   * first: a basis of the quotient ring as a vector space. Only for a
   * zero-dimensional ideal (otherwise there are infinitely many).
   */
  standardMonomials({ limit = 20000 } = {}) {
    if (this.isTrivial) return [];
    if (!this.isZeroDimensional) throw new RangeError("the quotient is not finite-dimensional: the ideal has dimension " + this.dimension);
    if (this._std) return this._std.map((e) => this._mono(e));
    const lms = this._basis.map((g) => g.lm);
    const n = this.ring.n;
    const seen = new Set();
    const out = [];
    const queue = [new Array(n).fill(0)];
    seen.add(keyOf(queue[0]));
    while (queue.length) {
      const e = queue.shift();
      out.push(e);
      if (out.length > limit) throw new RangeError(`more than ${limit} standard monomials`);
      for (let i = 0; i < n; i += 1) {
        const f = e.slice();
        f[i] += 1;
        const k = keyOf(f);
        if (seen.has(k) || lms.some((l) => dividesExp(l, f))) continue;
        seen.add(k);
        queue.push(f);
      }
    }
    out.sort(this._cmp);
    this._std = out;
    return out.map((e) => this._mono(e));
  }
  _mono(e) {
    return new MPoly(this.ring, new Map([[keyOf(e), { e, c: this.ring.K.one }]]));
  }
  /** The dimension of the quotient ring as a vector space: the number of solutions counted with multiplicity. */
  get quotientDimension() {
    if (this.isTrivial) return 0;
    if (!this.isZeroDimensional) return Infinity;
    this.standardMonomials();
    return this._std.length;
  }
  toString() {
    return `Gröbner basis (${this.ring.order}) of ${this.basis.length}: ${this.basis.join(", ")}`;
  }
}

/**
 * The reduced Gröbner basis of the ideal generated by some polynomials of
 * one ring, in that ring's monomial order (or { order } to use another).
 * { certificates: true } carries the cofactors from the start. maxPairs and
 * timeLimit (milliseconds) stop a computation that is too large, with a
 * RangeError.
 */
export function groebner(polys, { order, certificates = false, maxPairs = 200000, timeLimit = 0 } = {}) {
  const list = Array.isArray(polys) ? polys : [polys];
  if (!list.length) throw new RangeError("give at least one polynomial");
  const first = list.find((p) => p instanceof MPoly);
  if (!first) throw new TypeError("give polynomials of a ring: ring(\"x, y\").poly(\"x^2 - y\")");
  let R = first.ring;
  if (order !== undefined) R = R.withOrder(order);
  if (!(R instanceof PolyRing)) throw new TypeError("not a polynomial ring");
  void orderOf;
  const gens = list.map((p) => R.from(p)).filter((p) => !p.isZero());
  return new GroebnerBasis(R, gens, { certificates, maxPairs, timeLimit });
}
