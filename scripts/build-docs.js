// Build the docs site: bundle the library (with @ratmath/core) to
// docs/assets/structures.js, then write docs/<page>.html from docs/src/<page>.html.
//
// Each source page is a full HTML document. The build ends its <article> with
// the runnable code box (docs/assets/runbox.js, shared with the core, analytic
// and linalg docs; seed code from the page's <template data-runbox>), typesets
// \( … \) and \[ … \] with KaTeX, and stamps asset links with a content hash so
// a browser never mixes a cached script with a new page. Finally it runs every
// page's box against the library and fails if any of them throws.
//
//   bun run docs:build [--verbose]
import { readFileSync, writeFileSync, readdirSync } from "node:fs";
import { renderMath, seedOf, runboxHtml, stampAssets, stampLinks, KATEX_LINK } from "./docs-kit.js";
import { runCode } from "../docs/assets/runbox.js";
import { makeLang } from "../docs/assets/lang.js";
import { makeGroups } from "../docs/assets/groups.js";

const out = await Bun.build({ entrypoints: ["./scripts/docs-entry.js"], outdir: "./docs/assets", naming: "structures.js", target: "browser", format: "esm", minify: { whitespace: true, syntax: true, identifiers: true, keepNames: true } });
if (!out.success) {
  for (const m of out.logs) console.error(m);
  process.exit(1);
}
console.log("docs/assets/structures.js written");

const SITE_LINKS = `<footer class="sitelinks"><a href="https://github.com/ratmath/structures">GitHub</a> · <a href="https://ratmath.com/">ratmath.com</a> · <a href="https://ratmath.com/core/">core</a> · <a href="https://ratmath.com/shelf/">shelf</a> · <a href="https://ratmath.com/algebraic/">algebraic</a> · <a href="https://ratmath.com/linalg/">linalg</a> · <a href="https://ratmath.com/graphics/">graphics</a></footer>`;

const v = stampAssets("docs/assets");
const seeds = [];
const exprs = [];
const groupNames = [];
const attr = (s) => s.replace(/&quot;/g, '"').replace(/&lt;/g, "<").replace(/&gt;/g, ">").replace(/&amp;/g, "&");
let formulas = 0;
const files = readdirSync("docs/src").filter((f) => f.endsWith(".html")).sort();
for (const file of files) {
  const { seed, html } = seedOf(readFileSync(`docs/src/${file}`, "utf8"));
  seeds.push({ page: file, seed });
  // every expression a page offers (data-expr, data-examples, and chips listed in a data-check attribute)
  for (const m of html.matchAll(/data-(?:expr|examples|check)="([^"]*)"/g)) for (const x of attr(m[1]).split("|")) if (x.trim()) exprs.push({ page: file, expr: x });
  for (const m of html.matchAll(/data-groups?="([^"]*)"/g)) for (const x of attr(m[1]).split("|")) if (x.trim()) groupNames.push({ page: file, name: x });
  const box = runboxHtml(seed, { lib: "<code>@ratmath/structures</code> and Core's <code>Rational</code>" });
  if (!html.includes("</article>")) throw new Error(`${file}: no </article> to put the code box in`);
  const title = /<title>([^<]*)<\/title>/.exec(html)?.[1] ?? "RatMath Structures";
  const og = `<meta property="og:title" content="${title}">\n<meta property="og:description" content="Finite groups, polynomials in several variables, ideals, Gröbner bases and rings of numbers, exact, in JavaScript.">\n<meta property="og:type" content="website">\n<meta property="og:url" content="https://ratmath.com/structures/${file === "index.html" ? "" : file}">`;
  const math = renderMath(html.replace("</title>", `</title>\n${og}`).replace("</article>", `${box}\n${SITE_LINKS}\n</article>`));
  formulas += math.count;
  const page = math.count ? math.html.replace('<link rel="stylesheet" href="assets/style.css">', `${KATEX_LINK}\n<link rel="stylesheet" href="assets/style.css">`) : math.html;
  writeFileSync(`docs/${file}`, stampLinks(page, v));
}
console.log(`${files.length} pages written to docs/ (${formulas} formulas typeset, assets ?v=${v})`);

const scope = await import("./docs-entry.js");
let failed = 0;
for (const { page, seed } of seeds) {
  const res = await runCode(seed, scope);
  if (res.error) {
    failed++;
    console.error(`✗ ${page}: ${res.error.name}: ${res.error.message}`);
  } else if (process.argv.includes("--verbose")) {
    console.log(`— ${page} (${res.ms.toFixed(0)} ms)`);
    for (const pr of res.prints) console.log(`  ${pr.text}`);
    if (res.hasValue) console.log(`  → ${res.type ? `${res.type} ` : ""}${res.text}`);
  }
}
const lang = makeLang(scope);
for (const { page, expr } of exprs) {
  try {
    const env = {};
    for (const line of expr.split("\\n")) if (line.trim()) lang.run(line, env);
  } catch (e) {
    failed++;
    console.error(`✗ ${page}: ${expr}: ${e.message}`);
  }
}
const grp = makeGroups(scope);
for (const { page, name } of groupNames) {
  try {
    grp.report(grp.parse(name));
  } catch (e) {
    failed++;
    console.error(`✗ ${page}: group ${name}: ${e.message}`);
  }
}
console.log(`${exprs.length} page expressions and ${groupNames.length} group names checked`);
if (failed) {
  console.error(`${failed} example(s) failed`);
  process.exit(1);
}
console.log(`${seeds.length} code boxes run without errors`);
