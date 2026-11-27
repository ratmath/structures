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
        const args = name === "det" || name === "inverse" ? [] : splitTop(raw).map((a) => expand(a, env));
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
