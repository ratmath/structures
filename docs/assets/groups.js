// The group explorer used on the docs pages: reads a group's name such as
// "D4", "S4 × C2" or "GL(2,3)", and reports its invariants. No DOM here, so
// the docs build can check every name a page offers.

export function makeGroups(A) {
  const one = (s) => {
    const t = s.trim();
    let m;
    if ((m = /^\((.*)\)$/.exec(t))) return parse(m[1]);
    if ((m = /^C(\d+)$/.exec(t))) return A.cyclic(Number(m[1]));
    if ((m = /^D(\d+)$/.exec(t))) return A.dihedral(Number(m[1]));
    if ((m = /^Dic(\d+)$/.exec(t))) return A.dicyclic(Number(m[1]));
    if ((m = /^S(\d+)$/.exec(t))) return A.symmetric(Number(m[1]));
    if ((m = /^A(\d+)$/.exec(t))) return A.alternating(Number(m[1]));
    if (t === "Q8") return A.quaternion();
    if (t === "V4") return A.kleinFour();
    if ((m = /^Q(\d+)$/.exec(t)) && Number(m[1]) % 4 === 0) return A.dicyclic(Number(m[1]) / 4);
    if ((m = /^(GL|SL|PSL)\((\d+),\s*(\d+)\)$/.exec(t))) return A[m[1]](Number(m[2]), Number(m[3]));
    if ((m = /^Heis\((\d+)\)$/.exec(t))) return A.heisenberg(Number(m[1]));
    if ((m = /^SmallGroup\((\d+),\s*(\d+)\)$/.exec(t))) {
      const list = A.smallGroups(Number(m[1]));
      const G = list[Number(m[2]) - 1];
      if (!G) throw new Error(`there are ${list.length} groups of order ${m[1]}`);
      return G;
    }
    if ((m = /^C(\d+)\s*⋊\s*C(\d+)$/.exec(t))) throw new Error("write a semidirect product as SmallGroup(n, i), or build it in the code box with metacyclic(m, k, r)");
    throw new Error(`unknown group ${t}: try C6, D4, Dic3, Q8, V4, S4, A5, GL(2,3), SL(2,3), PSL(2,7), Heis(3), SmallGroup(16,4), or products like S3 × C2`);
  };
  function parse(text) {
    if (String(text).includes("|")) return A.presentation(String(text).trim()).toGroup({ maxCosets: 100000 });
    const parts = splitTop(String(text));
    const gs = parts.map(one);
    return gs.length === 1 ? gs[0] : A.directProduct(...gs);
  }

  /** The facts the explorer shows, as [label, text] rows, plus the group. */
  function report(G) {
    const rows = [];
    const yes = (b) => (b ? "yes" : "no");
    rows.push(["order", String(G.order)]);
    let id = null;
    try { id = G.order <= A.SMALL_ORDER_LIMIT ? A.identify(G) : null; } catch { id = null; }
    if (id) rows.push(["isomorphic to", `${id.name} (in the table of small groups)`]);
    rows.push(["generators", G.generators.map((g) => G.format(g)).join(", ") || "none"]);
    rows.push(["abelian / cyclic", `${yes(G.isAbelian)} / ${yes(G.isCyclic)}`]);
    rows.push(["element orders", [...G.orderStatistics()].map(([o, c]) => `${c} of order ${o}`).join(", ")]);
    rows.push(["class equation", String(G.classEquation())]);
    rows.push(["centre", `order ${G.centre().order}`]);
    rows.push(["commutator subgroup", `order ${G.derivedSubgroup().order}`]);
    rows.push(["solvable / nilpotent / simple", `${yes(G.isSolvable)} / ${yes(G.isNilpotent)} / ${yes(G.isSimple)}`]);
    rows.push(["normal subgroups (orders)", G.normalSubgroups().map((N) => N.order).join(", ")]);
    if (G.order <= 400) {
      const subs = G.subgroups();
      const by = new Map();
      for (const H of subs) by.set(H.order, (by.get(H.order) ?? 0) + 1);
      rows.push(["subgroups", `${subs.length}: ${[...by].map(([o, c]) => `${c} of order ${o}`).join(", ")}`]);
    }
    const ps = [...factor(G.order)];
    if (ps.length) rows.push(["Sylow counts", ps.map((p) => `n${p} = ${G.sylowCount(p)}`).join(", ")]);
    rows.push(["composition factors (orders)", G.compositionSeries().factors.join(", ") || "none"]);
    if (G.order <= 120) {
      const auts = A.automorphisms(G).length;
      rows.push(["automorphisms", `${auts} (inner: ${G.order / G.centre().order})`]);
    }
    return { group: G, rows };
  }

  return { parse, report };
}

/** Splits a data-groups list on "|", except inside ⟨ … ⟩ (a presentation). */
export function splitList(text) {
  const out = [];
  let depth = 0;
  let cur = "";
  for (const ch of String(text)) {
    if (ch === "⟨") depth += 1;
    if (ch === "⟩") depth -= 1;
    if (ch === "|" && depth === 0) { out.push(cur); cur = ""; continue; }
    cur += ch;
  }
  out.push(cur);
  return out.map((x) => x.trim()).filter(Boolean);
}

function factor(n) {
  const out = new Set();
  for (let p = 2; p * p <= n; p += 1) while (n % p === 0) { out.add(p); n /= p; }
  if (n > 1) out.add(n);
  return out;
}

function splitTop(s) {
  const out = [];
  let depth = 0;
  let cur = "";
  for (const ch of s.replace(/\s+x\s+/g, " × ")) {
    if (ch === "(") depth += 1;
    if (ch === ")") depth -= 1;
    if ((ch === "×" || ch === "*") && depth === 0) {
      out.push(cur);
      cur = "";
      continue;
    }
    cur += ch;
  }
  out.push(cur);
  return out.filter((p) => p.trim());
}
