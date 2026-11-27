// The small expression language of the docs pages and the explorer.
//
//   x^2 y - 3/2 x z + 1        polynomials and fractions, written as on paper
//   f = (x + y)^3              keep a value under a name
//   factor(f), gcd(f, g), resultant(f, g, y), det(a, b; c, d) …
//   order lex                  change the monomial order
//
// Variables are single letters; a ring is made from whichever letters occur.
// makeLang(A) takes the library (so the build can check every example).

export function makeLang(A) {
  const letters = (s) => [...new Set(s.match(/[A-Za-z]/g) ?? [])].sort();

  /** Split at commas (or another separator) that are not inside parentheses. */
  function splitTop(s, sep = ",") {
    const out = [];
    let depth = 0;
    let cur = "";
    for (const ch of s) {
      if (ch === "(") depth += 1;
      else if (ch === ")") depth -= 1;
      if (ch === sep && depth === 0) {
        out.push(cur);
        cur = "";
      } else cur += ch;
    }
    out.push(cur);
    return out.map((x) => x.trim());
  }

  const text = (t) => ({ type: "text", text: t });

  /** A context: one ring holding every letter of the given expanded texts. */
  function context(env, texts, extra = []) {
    const vars = [...new Set([...texts.flatMap(letters), ...extra])].sort();
    const R = A.ring(vars.length ? vars : ["x"], { order: env.order ?? "grevlex" });
    const F = A.fractions(R);
    const frac = (t) => F.parse(t);
    const poly = (t) => {
      const v = F.parse(t);
      if (!v.isPolynomial()) throw new Error(`${v} is not a polynomial`);
      return v.toPolynomial();
    };
    return { R, F, frac, poly };
  }

  const variable = (s) => {
    if (!/^[A-Za-z]$/.test(s)) throw new Error(`expected a variable (one letter), got "${s}"`);
    return s;
  };
  const whole = (s) => {
    if (!/^\d+$/.test(s)) throw new Error(`expected a whole number, got "${s}"`);
    return Number(s);
  };
  const int = (s) => {
    const t = s.replace(/[\s()]/g, "").replace("−", "-");
    if (!/^-?\d+$/.test(t)) throw new Error(`expected a whole number, got "${s}"`);
    return Number(t);
  };
  const factored = (unit, factors) => {
    const u = String(unit);
    const parts = factors.map(({ factor, multiplicity }) => {
      const s = String(factor);
      const body = factor.length > 1 || factors.length > 1 || u !== "1" ? `(${s})` : s;
      return multiplicity > 1 ? `${body}${String(multiplicity).replace(/\d/g, (d) => "⁰¹²³⁴⁵⁶⁷⁸⁹"[d])}` : body;
    });
    if (u !== "1" || !parts.length) parts.unshift(u === "-1" && parts.length ? "−" : u.replace("-", "−"));
    return parts.join(" · ").replace(/^− · /, "−");
  };
  const matrixOf = (raw, env) => splitTop(raw, ";").map((row) => splitTop(row).map((x) => expand(x, env)));

  function inverse(K, m) {
    const n = m.length;
    const a = m.map((row, i) => [...row, ...Array.from({ length: n }, (_, j) => (i === j ? K.one : K.zero))]);
    for (let c = 0; c < n; c += 1) {
      const p = a.findIndex((row, i) => i >= c && K.sign(row[c]) !== 0);
      if (p < 0) throw new Error("this matrix has no inverse: its determinant is 0");
      [a[c], a[p]] = [a[p], a[c]];
      const lead = a[c][c];
      a[c] = a[c].map((v) => K.div(v, lead));
      for (let i = 0; i < n; i += 1) {
        if (i === c) continue;
        const k = a[i][c];
        if (K.sign(k) !== 0) a[i] = a[i].map((v, j) => K.sub(v, K.mul(k, a[c][j])));
      }
    }
    return a.map((row) => row.slice(n));
  }

  const real = (z) => (typeof z === "number" ? z : z.re);
  const sub = (n) => String(n).replace(/\d/g, (d) => "₀₁₂₃₄₅₆₇₈₉"[d]);
  /** tower(a^2 - 2, b^3 - a): each polynomial introduces its newest letter. */
  function towerOf(raw) {
    const defs = {};
    for (const part of splitTop(raw)) {
      const fresh = letters(part).filter((ch) => !Object.hasOwn(defs, ch));
      if (fresh.length !== 1) throw new Error(`"${part}" should bring in exactly one new letter`);
      defs[fresh[0]] = part;
    }
    return A.extension(defs);
  }
  const RAW = new Set(["det", "inverse", "large", "tower", "galois", "radicals", "unity", "polygon", "constructible", "conway"]);

  // Each function gets its arguments already expanded to plain text (names
  // replaced), plus the raw inside of the parentheses for matrices. It returns
  // text that can be read again (a value) or a result to show.
  const FUNCS = {
    expand: ([f], env) => context(env, [f]).frac(f).toText(),
    factor: ([f], env) => {
      const { unit, factors } = A.factor(context(env, [f]).poly(f));
      return text(factored(unit, factors));
    },
    sqfree: ([f], env) => {
      const { unit, factors } = A.squareFree(context(env, [f]).poly(f));
      return text(factored(unit, factors));
    },
    irreducible: ([f], env) => text(A.isIrreducible(context(env, [f]).poly(f)) ? "irreducible over Q" : "not irreducible over Q"),
    gcd: (args, env) => A.gcd(args.map(context(env, args).poly)).toText(),
    lcm: (args, env) => A.lcm(args.map(context(env, args).poly)).toText(),
    content: ([f, v], env) => A.content(context(env, [f], [variable(v)]).poly(f), v).toText(),
    primpart: ([f, v], env) => A.primitivePart(context(env, [f], [variable(v)]).poly(f), v).toText(),
    diff: ([f, v, n = "1"], env) => {
      let r = context(env, [f], [variable(v)]).frac(f);
      for (let k = 0; k < whole(n); k += 1) r = r.diff(v);
      return r.toText();
    },
    subs: ([f, v, x], env) => {
      const c = context(env, [f, x], [variable(v)]);
      return c.frac(f).subs({ [v]: c.poly(x) }).toText();
    },
    resultant: ([f, g, v], env) => {
      const c = context(env, [f, g], [variable(v)]);
      return A.resultant(c.poly(f), c.poly(g), v).toText();
    },
    disc: ([f, v], env) => A.discriminant(context(env, [f], [variable(v)]).poly(f), v).toText(),
    degree: ([f, v], env) => {
      const p = context(env, [f], v ? [variable(v)] : []).poly(f);
      return String(v ? p.degree(v) : p.degree());
    },
    lead: ([f], env) => context(env, [f]).poly(f).leadingTerm().toText(),
    homogenize: ([f], env) => context(env, [f]).poly(f).homogenize("h").toText(),
    rem: ([f, ...gs], env) => {
      const c = context(env, [f, ...gs]);
      return c.poly(f).rem(gs.map(c.poly)).toText();
    },
    divide: ([f, ...gs], env) => {
      const c = context(env, [f, ...gs]);
      const divisors = gs.map(c.poly);
      const { quotients, remainder } = c.poly(f).divmod(divisors);
      return {
        type: "record",
        title: `remainder ${remainder}`,
        rows: [...divisors.map((g, i) => [`times ${g}`, String(quotients[i])]), ["order", c.R.order]],
      };
    },
    det: (_, env, raw) => {
      const m = matrixOf(raw, env);
      const c = context(env, m.flat());
      return A.det(m.map((row) => row.map(c.poly)), c.R).toText();
    },
    inverse: (_, env, raw) => {
      const m = matrixOf(raw, env);
      if (m.some((row) => row.length !== m.length)) throw new Error("the inverse needs a square matrix: rows are separated by ;");
      const c = context(env, m.flat());
      return { type: "matrix", rows: inverse(c.F, m.map((row) => row.map(c.frac))).map((row) => row.map(String)) };
    },
    jacobian: (args, env) => {
      const c = context(env, args);
      return { type: "matrix", head: c.R.vars.map((v) => `∂/∂${v}`), rows: A.jacobian(args.map(c.poly)).map((row) => row.map(String)) };
    },
    hessian: ([f], env) => ({ type: "matrix", rows: A.hessian(context(env, [f]).poly(f)).map((row) => row.map(String)) }),
    elementary: ([k, ...vs]) => A.elementary(A.ring(vs.map(variable)), whole(k)).toText(),
    powersum: ([k, ...vs]) => A.powerSum(A.ring(vs.map(variable)), whole(k)).toText(),
    complete: ([k, ...vs]) => A.completeHomogeneous(A.ring(vs.map(variable)), whole(k)).toText(),
    groebner: (args, env) => {
      const c = context(env, args);
      const G = A.groebner(args.map(c.poly));
      return { type: "record", title: G.isTrivial ? "the whole ring: no common zeros" : `${G.length} polynomial${G.length === 1 ? "" : "s"} (${c.R.order}, ${c.R.vars.join(" > ")})`, rows: [...G.basis.map((g, i) => [`g${i + 1}`, String(g)]), ["solution set", G.dimension < 0 ? "empty" : G.dimension === 0 ? `${G.quotientDimension} point${G.quotientDimension === 1 ? "" : "s"}, counted with multiplicity` : `dimension ${G.dimension}`]] };
    },
    member: ([f, ...gs], env) => {
      const c = context(env, [f, ...gs]);
      const I = A.ideal(c.R, gs.map(c.poly));
      const F = c.poly(f);
      const r = I.reduce(F);
      if (!r.isZero()) return { type: "record", title: "not in the ideal", rows: [["remainder", String(r)], ["some power in it?", I.radicalContains(F) ? `yes: the power ${I.exponentOf(F)}` : "no"]] };
      const ex = I.express(F);
      return { type: "record", title: "in the ideal", rows: I.generators.map((g, i) => [`times ${g}`, String(ex.cofactors[i])]) };
    },
    eliminate: ([v, ...gs], env) => {
      const vs = [...v.replace(/\s/g, "")].map(variable);
      const c = context(env, gs, vs);
      const J = A.ideal(c.R, gs.map(c.poly)).eliminate(vs);
      return { type: "record", title: J.isZero ? "nothing: no relation is left" : `without ${vs.join(", ")}`, rows: J.basis.map((g, i) => [`h${i + 1}`, String(g)]) };
    },
    radical: (args, env) => {
      const c = context(env, args);
      const I = A.ideal(c.R, args.map(c.poly));
      const J = I.radical();
      return { type: "record", title: I.equals(J) ? "already radical" : "the radical is larger", rows: J.basis.map((g, i) => [`r${i + 1}`, String(g)]) };
    },
    solve: (args, env) => {
      const c = context(env, args);
      const s = A.solve(c.R, args.map(c.poly), { complex: true });
      if (!s.finite) return text(`infinitely many solutions: the solution set has dimension ${s.dimension}`);
      if (!s.count) return text("no solutions, even over the complex numbers");
      const show = (p) => c.R.vars.map((v) => `${v} = ${p.values[v]}`).join(",  ");
      const rows = [];
      s.real.forEach((p, i) => rows.push([`real ${i + 1}`, `${show(p)}   ≈ (${c.R.vars.map((v) => p.approx[v].toFixed(6)).join(", ")})`]));
      s.complex.slice(0, 12).forEach((p, i) => rows.push([`complex ${i + 1}`, show(p)]));
      if (s.complex.length > 12) rows.push(["", `and ${s.complex.length - 12} more complex solutions`]);
      rows.push(["with multiplicity", String(s.multiplicity)]);
      return { type: "record", title: `${s.count} solution${s.count === 1 ? "" : "s"}, ${s.real.length} real`, rows };
    },
    zmod: ([n]) => {
      const Z = A.Zn(whole(n));
      if (Z.n > 5000n) throw new Error("zmod lists elements, so keep n up to 5000");
      const list = (xs) => (xs.length ? (xs.length > 40 ? `${xs.slice(0, 40).join(", ")}, … (${xs.length})` : xs.join(", ")) : "none");
      const d = Z.decompose();
      return { type: "record", title: `${Z}: ${Z.isField ? "a field" : Z.isReduced ? "not a domain, no nilpotents" : "not a domain"}`, rows: [["units", `${list(Z.units())}`], ["unit group", `${Z.unitGroup().text}, order ${Z.unitCount}${Z.primitiveRoot() === null ? "" : `, generated by ${Z.primitiveRoot()}`}`], ["zero divisors", list(Z.zeroDivisors())], ["nilpotents", list(Z.nilpotents())], ["idempotents", list(Z.idempotents())], ["ideals", Z.ideals().map((I) => `(${I.generator})`).join(", ")], ["splits as", d.moduli.map((m) => `Z/${m}`).join(" × ")]] };
    },
    crt: (args) => {
      if (!args.length || args.length % 2) throw new Error("crt takes pairs: crt(2, 3, 3, 5) solves x ≡ 2 mod 3 and x ≡ 3 mod 5");
      const rs = args.filter((_, i) => i % 2 === 0).map(int);
      const ms = args.filter((_, i) => i % 2 === 1).map(whole);
      const r = A.crt(rs, ms);
      return text(r ? `x ≡ ${r.value} mod ${r.modulus}` : "no solution: the congruences contradict each other");
    },
    quadratic: ([d]) => {
      const K = A.quadraticRing(int(d));
      const cg = K.classGroup();
      const u = K.units;
      const small = [2, 3, 5, 7, 11, 13].map((p) => `${p} ${K.splitting(p).kind}`).join(", ");
      return { type: "record", title: `${K}`, rows: [["discriminant", String(K.discriminant)], ["units", u.finite ? u.list.join(", ") : `${u.text}, norm ${u.norm}`], ["class number", `${cg.classNumber}${cg.classNumber > 1 ? ` (${cg.text})` : ""}`], ["unique factorisation", K.isUFD ? "yes" : "no"], ["division with remainder by the norm", K.isNormEuclidean ? "yes" : "no"], ["ideal classes", cg.representatives.map(String).join(", ")], ["small primes", small]] };
    },
    qfactor: ([d, a, b = "0"]) => {
      const K = A.quadraticRing(int(d));
      const x = K.sqrt(int(a), int(b));
      const show = (f) => `${String(f.unit) === "1" ? "" : `(${f.unit}) · `}${f.factors.map((q) => `(${q})`).join(" · ")}`;
      const fs = K.factorizations(x);
      const rows = fs.map((f, i) => [fs.length > 1 ? `way ${i + 1}` : "irreducibles", show(f)]);
      rows.push(["as ideals", K.ideal(x).factor().map(({ prime, exponent }) => `${prime}${exponent > 1 ? `^${exponent}` : ""}`).join(" · ") || "(1)"]);
      return { type: "record", title: `${x} in ${K}: ${fs.length === 1 ? "one factorisation" : `${fs.length} different factorisations`}`, rows };
    },
    large: (_, env, raw) => {
      const name = raw.trim();
      const L = A.largeGroups;
      const m = /^([SAM])\s*(\d+)$/i.exec(name);
      let G;
      if (/^rubik/i.test(name)) G = L.rubik();
      else if (/^pocket/i.test(name)) G = L.pocketCube();
      else if (m) {
        const n = Number(m[2]);
        const kind = m[1].toUpperCase();
        if (kind !== "M" && n > 60) throw new Error("degrees up to 60 here, please");
        G = kind === "S" ? L.symmetric(n) : kind === "A" ? L.alternating(n) : L.mathieu(n);
      } else throw new Error(`large(…) takes rubik, pocket, S n, A n or M11, M12, M22, M23, M24, not "${name}"`);
      const orbits = G.orbits().filter((o) => o.length > 1).map((o) => o.length);
      const rows = [
        ["acts on", `${G.degree} points`],
        ["order", String(G.order)],
        ["as a product", G.orbitLengths.join(" · ")],
        ["base", G.base.join(", ")],
        ["strong generators", String(G.strongGenerators.length)],
        ["orbits", orbits.join(" + ") || "none"],
      ];
      if (G.isTransitive) rows.push(["primitive", G.isPrimitive ? "yes" : `no: blocks of size ${G.minimalBlocks()[0].length}`]);
      if (G.degree <= 24) {
        const D = G.derivedSubgroup();
        rows.push(["commutator subgroup", `order ${D.order}${D.order === G.order ? " (perfect)" : ""}`], ["solvable", G.isSolvable ? "yes" : "no"]);
      }
      return { type: "record", title: String(G), rows };
    },
    tower: (_, env, raw) => {
      const K = towerOf(raw);
      const prim = K.primitiveElement();
      const rows = [
        ["degree", K.towerLaw()],
        ["basis over Q", K.basis.map(String).join(", ")],
        ["primitive element", `${prim.element}, a root of ${prim.text}`],
      ];
      if (K.degree <= 8) {
        const auts = K.automorphisms();
        rows.push(["automorphisms", `${auts.length}${auts.length === K.degree ? ", so it is a Galois extension" : `, fewer than the degree ${K.degree}: not Galois`}`]);
        for (const s of auts.slice(1, 9)) rows.push(["", String(s)]);
      }
      return { type: "record", title: K.describe(), rows };
    },
    galois: (_, env, raw) => {
      const K = towerOf(raw);
      if (K.degree > 8) throw new Error("degree up to 8 here, please");
      if (!K.isGalois) return text(`${K} has ${K.automorphisms().length} automorphisms but degree ${K.degree}, so it is not a Galois extension of Q`);
      const G = K.galoisGroup();
      const rows = [["Galois group", `${A.identify(G).name}, of order ${G.order}`]];
      const subs = [...G.subgroups()].sort((x, y) => x.order - y.order);
      for (const H of subs) {
        const F = K.fixedField([...H.elements]);
        rows.push([`subgroup of order ${H.order}${G.isNormal(H) ? " (normal)" : ""}`, `fixes a field of degree ${F.degree}: ${F.degree === 1 ? "Q" : `Q(${F.generator}), ${F.text} = 0`}`]);
      }
      return { type: "record", title: `${K}: subgroups and the fields they fix`, rows };
    },
    radicals: (_, env, raw) => {
      const r = A.solveByRadicals(raw.trim());
      const rows = [];
      for (const f of r.factors) {
        rows.push([r.factors.length > 1 ? f.text : "Galois group", r.factors.length > 1 ? `Galois group ${f.group}, of order ${f.order}` : `${f.group}, of order ${f.order}`]);
        if (!f.solvable) {
          rows.push(["", "not solvable: no formula in radicals"]);
          continue;
        }
        if (f.chain.length) rows.push(["steps", f.chain.join(", ")]);
        const shown = new Set();
        f.roots.forEach((root, i) => {
          const e = root.expression;
          const line = e.inline({ maxLength: 160 });
          if (line !== null) rows.push([`x${sub(i + 1)}`, line]);
          else {
            rows.push([`x${sub(i + 1)}`, e.text]);
            for (const d of e.definitions) if (!shown.has(d.text)) (shown.add(d.text), rows.push(["", d.text]));
          }
        });
        const ok = f.roots.every((root) => root.expression.check(30).agrees);
        rows.push(["check", ok ? "every formula agrees with its root to 30 digits" : "a formula failed its check"]);
      }
      return { type: "record", title: r.solvable ? `${r.text} is solvable by radicals` : `${r.text} is not solvable by radicals`, rows };
    },
    unity: (_, env, raw) => {
      const n = whole(raw.trim());
      if (n > 60) throw new Error("n up to 60 here, please");
      const e = A.rootOfUnityByRadicals(n);
      const z = e.approx();
      const rows = [];
      const line = e.inline({ maxLength: 200 });
      if (line !== null) rows.push([`ζ${sub(n)}`, line]);
      else e.definitions.forEach((d, i) => rows.push([i ? "" : "step by step", d.text]));
      rows.push(["value", `${z.re.toFixed(12)} ${z.im < 0 ? "−" : "+"} ${Math.abs(z.im).toFixed(12)}·i`], ["radicals used", [...new Set(e.indices)].sort((a, b) => a - b).map((k) => (k === 2 ? "square roots" : k === 3 ? "cube roots" : `${k}th roots`)).join(", ") || "none"], ["check", e.check(30).agrees ? "agrees with e^(2πi/n) to 30 digits" : "failed"]);
      return { type: "record", title: `e^(2πi/${n}) in radicals`, rows };
    },
    polygon: (_, env, raw) => {
      const n = whole(raw.trim());
      if (n > 300) throw new Error("n up to 300 here, please");
      const p = A.regularPolygon(n);
      const rows = [["φ(n)", String(p.totient)], ["why", p.reason]];
      if (p.constructible) {
        const line = p.cos.inline({ maxLength: 220 });
        if (line !== null) rows.push([`cos(2π/${n})`, line]);
        else {
          rows.push([`cos(2π/${n})`, p.cos.text]);
          const defs = p.cos.definitions;
          for (const d of defs.slice(0, 12)) rows.push(["", d.text]);
          if (defs.length > 12) rows.push(["", `… and ${defs.length - 12} more square roots`]);
        }
        rows.push(["value", real(p.cos.approx()).toFixed(15)]);
      }
      return { type: "record", title: `The regular ${n}-gon ${p.constructible ? "can" : "cannot"} be drawn with ruler and compass`, rows };
    },
    constructible: (_, env, raw) => {
      const r = A.constructible(raw.trim());
      const rows = r.factors.map((f) => [f.text, `${f.group ? `Galois group ${f.group}, of order ${f.order}` : `degree ${f.degree}`}: ${f.constructible ? "constructible" : "not constructible"}`]);
      rows.push(["why", r.reason]);
      const word = r.constructible === true ? "Every root can" : r.constructible === "some" ? "Some roots can" : r.constructible === false ? "No root can" : "Undecided whether the roots can";
      return { type: "record", title: `${word} be constructed with ruler and compass`, rows };
    },
    conway: (_, env, raw) => {
      const [p, n] = splitTop(raw).map(whole);
      if (n < 1 || p ** n > 5000) throw new Error("fields up to 5000 elements here, please");
      const F = A.conwayField(p, n);
      const rows = [["Conway polynomial", F.text], ["field", `${F.size} elements, generator a of order ${BigInt(F.size) - 1n}`]];
      for (let m = 1; m < n; m += 1) {
        if (n % m) continue;
        const e = F.embed(m);
        rows.push([`GF(${p}${m > 1 ? `^${m}` : ""}) inside`, `its generator is a^${e.exponent} = ${F.format(e.image)}`]);
      }
      return { type: "record", title: `GF(${p}^${n}) = GF(${p})[a] / (${F.text.replaceAll("x", "a")})`, rows };
    },
    symmetric: ([f, ...vs], env) => {
      const c = context(env, [f], vs.map(variable));
      const p = c.poly(f);
      if (!A.isSymmetric(p)) return text(`not symmetric in ${c.R.vars.join(", ")}`);
      const { polynomial } = A.toElementary(p);
      return { type: "record", title: String(polynomial), rows: c.R.vars.map((_, i) => [`e${i + 1}`, String(A.elementary(c.R, i + 1))]) };
    },
  };

  /** Replace names and function calls, leaving text the polynomial reader takes. */
  function expand(src, env, top = false) {
    let out = "";
    let i = 0;
    while (i < src.length) {
      const m = /^[A-Za-z][A-Za-z0-9]*/.exec(src.slice(i));
      if (!m) {
        out += src[i];
        i += 1;
        continue;
      }
      let name = m[0];
      if (!Object.hasOwn(FUNCS, name) && !Object.hasOwn(env.names, name)) name = /^[A-Za-z]+/.exec(name)[0];
      const after = src.slice(i + name.length);
      const open = /^\s*\(/.exec(after);
      if (Object.hasOwn(FUNCS, name) && open) {
        const start = i + name.length + open[0].length;
        let depth = 1;
        let j = start;
        for (; j < src.length && depth > 0; j += 1) {
          if (src[j] === "(") depth += 1;
          else if (src[j] === ")") depth -= 1;
        }
        if (depth > 0) throw new Error(`${name}( is never closed`);
        const raw = src.slice(start, j - 1);
        const args = RAW.has(name) ? [] : splitTop(raw).map((a) => expand(a, env));
        const value = FUNCS[name](args, env, raw);
        const alone = src.slice(0, i).trim() === "" && src.slice(j).trim() === "";
        if (typeof value !== "string") {
          if (top && alone) return value;
          throw new Error(`${name}(…) gives a table, so it has to stand on its own`);
        }
        out += `(${value})`;
        i = j;
        continue;
      }
      if (Object.hasOwn(env.names, name)) out += `(${env.names[name]})`;
      else out += [...name].map((ch) => (Object.hasOwn(env.names, ch) ? `(${env.names[ch]})` : ch)).join(" ");
      i += name.length;
    }
    return out;
  }

  /**
   * Run one line. `env` keeps names and the order between lines. Returns
   * { name, result } with result one of { type: "text" | "record" | "matrix" }.
   */
  function run(line, env = {}) {
    env.names ??= {};
    const src = line.trim();
    const ord = /^order\s*=?\s*(\w+)$/.exec(src);
    if (ord) {
      if (!Object.hasOwn(A.ORDERS, ord[1])) throw new Error(`unknown order "${ord[1]}": use lex, grlex or grevlex`);
      env.order = ord[1];
      return { name: null, result: text(`terms are now ordered by ${ord[1]}`) };
    }
    const assign = /^([A-Za-z][A-Za-z0-9]*)\s*=(?!=)\s*(.+)$/s.exec(src);
    if (assign && !Object.hasOwn(FUNCS, assign[1])) {
      const value = expand(assign[2], env, true);
      if (typeof value !== "string") throw new Error("only polynomials and fractions can be kept under a name");
      const v = context(env, [value]).frac(value);
      env.names[assign[1]] = v.toText();
      return { name: assign[1], result: text(String(v)) };
    }
    const value = expand(src, env, true);
    if (typeof value !== "string") return { name: null, result: value };
    return { name: null, result: text(String(context(env, [value]).frac(value))) };
  }

  return { run, functions: Object.keys(FUNCS) };
}
