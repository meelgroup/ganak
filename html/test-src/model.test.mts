// Tests of the "Probability models" tab (model.js) against ganak-model.js,
// ganak.js and bdnnf.js. `npm test` in html/.

import { test } from "node:test";
import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import { spawnSync } from "node:child_process";
import path from "node:path";
import {
  MODEL_EXAMPLES, Q, fmtP, parseLits, batchFile, answers, impacts,
  modelDot, pretty, gmTokens, gmHighlight, type GmModel,
} from "../model.js";
import { HTML, runTool } from "./tools.mjs";

const dec = (b: Uint8Array) => new TextDecoder().decode(b);

async function translate(text: string) {
  const r = await runTool("ganak-model", ["/m.gm", "-o", "/m.cnf", "--json", "/m.json"], { "/m.gm": text },
    ["/m.cnf", "/m.json"]);
  const model: GmModel = JSON.parse(dec(r.files["/m.json"]));
  return { r, model, cnf: r.files["/m.cnf"] ? dec(r.files["/m.cnf"]) : "" };
}

async function compiled(text: string) {
  const t = await translate(text);
  assert.equal(t.r.status, 0, t.r.stderr);
  const c = await runTool("ganak", ["--verb", "0", "--mode", "1", "--compile", "/out.bdnnf", "/m.cnf"],
    { "/m.cnf": t.cnf }, ["/out.bdnnf"]);
  assert.equal(c.status, 0, c.stderr);
  const bd = (args: string[], extra: Record<string, string> = {}) =>
    runTool("bdnnf", args, { "/c.bdnnf": c.files["/out.bdnnf"], "/m.cnf": t.cnf, ...extra });
  return { ...t, bd };
}

const HYD = MODEL_EXAMPLES.find((e) => e.name.startsWith("Airliner"))!.text;

test("Q: exact arithmetic and display", () => {
  const a = Q.parse("1/3"), b = Q.parse("2/6");
  assert.equal(a.cmp(b), 0);
  assert.equal(a.add(Q.parse("1/6")).toString(), "1/2");
  assert.equal(Q.one.sub(Q.parse("1/50000")).toString(), "49999/50000");
  const tiny = Q.parse("1260800407259559027271693566674542840180999/500000000000000000000000000000000000000000000000000");
  assert.ok(Math.abs(tiny.toNumber() / 2.5216008e-9 - 1) < 1e-6);
  assert.equal(fmtP(tiny), "2.52×10⁻⁹");
  assert.equal(fmtP(Q.parse("287/290")), "0.9897");
  assert.equal(fmtP(null), "undefined");
  assert.equal(fmtP(Q.one.sub(tiny)), "1 − 2.52×10⁻⁹");
  assert.equal(fmtP(Q.parse("9995/10000")), "1 − 5.00×10⁻⁴");
});

test("the hydraulics example is examples/models/hydraulics.gm", () => {
  assert.equal(HYD, readFileSync(path.join(HTML, "..", "examples", "models", "hydraulics.gm"), "utf8"));
});

for (const ex of MODEL_EXAMPLES)
  test(`example "${ex.name}": the page's answers equal ganak-model --solve`, async () => {
    const c = await compiled(ex.text);
    assert.deepEqual(c.model.errors, []);
    const b = await c.bd(["batch", "/c.bdnnf", "/q.txt", "--cnf", "/m.cnf"], { "/q.txt": batchFile(c.model) });
    assert.equal(b.status, 0, b.stderr);
    const page = answers(c.model, b.stdout);
    const s = await runTool("ganak-model", ["--solve", "/m.gm"], { "/m.gm": ex.text });
    const lines = s.stdout.split("\n").filter((l) => l.startsWith("P(") && !l.startsWith("P(constraints)") || l.startsWith("require"));
    assert.equal(lines.length, page.length);
    page.forEach((a, i) => {
      const exact = lines[i].match(/= (\S+) ~/) ?? lines[i].match(/\((\S+) ~/);
      assert.ok(exact, lines[i]);
      assert.equal(a.p!.toString(), exact![1], pretty(a.stmt));
      if (a.stmt.kind === "require") assert.equal(a.pass, lines[i].includes("PASS"));
    });
  });

test("hydraulics: answers, marginals and what matters", async () => {
  const c = await compiled(HYD);
  const b = await c.bd(["batch", "/c.bdnnf", "/q.txt", "--cnf", "/m.cnf"], { "/q.txt": batchFile(c.model) });
  const ans = answers(c.model, b.stdout);
  assert.equal(pretty(ans[0].stmt), "P(lost)");
  assert.equal(fmtP(ans[0].p), "2.52×10⁻⁹");
  assert.equal(ans[ans.length - 1].pass, false);

  const node = (n: string) => c.model.nodes.find((x) => x.name === n)!;
  const lost = node("lost"), leak = node("leak_g");
  const base = parseLits((await c.bd(["lits", "/c.bdnnf", "--cnf", "/m.cnf"])).stdout);
  const withTop = parseLits((await c.bd(["lits", "/c.bdnnf", "--cnf", "/m.cnf", "--assume", String(lost.var)])).stdout);
  const { pT, rows } = impacts(c.model, base, withTop);
  assert.equal(pT!.toString(), ans[0].p!.toString());
  const count = async (lits: number[]) =>
    Q.parse((await c.bd(["count", "/c.bdnnf", "--cnf", "/m.cnf", "--assume", lits.join(",")])).stdout.trim());
  const x = rows.find((r) => r.node.name === "leak_g")!;
  assert.equal(x.ifYes!.toString(), (await count([leak.var, lost.var])).div(await count([leak.var])).toString());
  assert.equal(x.ifNo!.toString(), (await count([-leak.var, lost.var])).div(await count([-leak.var])).toString());
  assert.equal(x.p.toString(), "1/20000");
  assert.equal(rows.length, 13);
  const swing = (r: typeof x) => Math.abs(r.ifYes!.toNumber() - r.ifNo!.toNumber());
  for (let i = 1; i < rows.length; i++) assert.ok(swing(rows[i]) <= swing(rows[i - 1]) * (1 + 1e-9));
  assert.equal(rows[0].node.name, "dual");

  const marg = new Map<string, Q>();
  const z = (base.get(1) ?? Q.zero).add(base.get(-1) ?? Q.zero);
  for (const n of c.model.nodes) marg.set(n.name, base.get(n.var)!.div(z));
  const dot = modelDot(c.model, marg, marg, new Map([["eng1", true]]));
  assert.match(dot, /n\d+ \[id="gm\d+", shape=box, style="filled,rounded"/);
  const lostLine = dot.split("\n").find((l) => l.includes("Loss of control"))!;
  assert.match(lostLine, /shape=box, style="filled", peripheries=2/);
  const ptuLine = dot.split("\n").find((l) => l.includes("Power transfer unit fails"))!;
  assert.match(ptuLine, /shape=ellipse/);
  assert.match(dot, /penwidth=3, color="#16a34a"/);
  assert.ok(dot.includes("Loss of control"));
});

test("constraints are drawn, linked to what they mention", async () => {
  const t = await translate(MODEL_EXAMPLES.find((e) => e.name.startsWith("Weather"))!.text);
  const c = t.model.stmts.find((s) => s.kind === "constraint")!;
  assert.deepEqual(c.deps!.map((d) => d.name).sort(), ["rain", "sun"]);
  const dot = modelDot(t.model, null, null, new Map());
  assert.match(dot, /c\d+ \[id="gc\d+", shape=hexagon.*label=<!\(rain &amp; sun\)>/);
  assert.equal((dot.match(/shape=hexagon/g) ?? []).length, 4);
});

// The page renders with viz.js, which is Graphviz: it rejects the same labels as dot.
const HAS_DOT = spawnSync("dot", ["-V"]).status === 0;

function renderSvg(dot: string): string {
  const r = spawnSync("dot", ["-Tsvg"], { input: dot, encoding: "utf8" });
  assert.equal(r.status, 0, r.stderr);
  assert.equal(r.stderr, "", r.stderr);
  return r.stdout;
}

function assertNodesLabelled(svg: string, m: GmModel) {
  const groups = svg.split('<g id="gm').slice(1);
  assert.equal(groups.length, m.nodes.length);
  for (const g of groups) assert.ok(g.slice(0, g.indexOf("</g>")).includes("<text"), "unlabelled node gm" + g.slice(0, 4));
}

test("impossible evidence: every node keeps its name, no probabilities", async () => {
  const c = await compiled(HYD);
  const v = (n: string) => c.model.nodes.find((x) => x.name === n)!.var;
  const r = await c.bd(["count", "/c.bdnnf", "--cnf", "/m.cnf", "--assume", `${v("lost")},${-v("land_asap")}`]);
  assert.equal(Q.parse(r.stdout.trim()).isZero(), true);
  const dot = modelDot(c.model, null, null, new Map([["lost", true], ["land_asap", false]]), true);
  assert.ok(!dot.includes("<FONT POINT-SIZE=\"10\"></FONT>"));
  const lines = dot.split("\n").filter((l) => /^\s+n\d+ \[id=/.test(l));
  assert.equal(lines.length, c.model.nodes.length);
  for (const l of lines) {
    assert.ok(l.includes("P undefined"), l);
    assert.ok(!/">[pP] = /.test(l), l);
    assert.ok(l.includes('fillcolor="#f1f3f5"'), l);
  }
  assert.match(lines.find((l) => l.includes("Loss of control"))!, /penwidth=3, color="#16a34a"/);
  assert.match(lines.find((l) => l.includes("Land at the nearest airport"))!, /penwidth=3, color="#dc2626"/);
  if (HAS_DOT) {
    const svg = renderSvg(dot);
    assertNodesLabelled(svg, c.model);
    assert.equal((svg.match(/<g id="gm\d+" class="node ev">/g) ?? []).length, 2);
  }
});

for (const ex of MODEL_EXAMPLES)
  test(`example "${ex.name}": the graph renders, every node labelled`, { skip: !HAS_DOT && "no dot" }, async () => {
    const c = await compiled(ex.text);
    const base = parseLits((await c.bd(["lits", "/c.bdnnf", "--cnf", "/m.cnf"])).stdout);
    const z = (base.get(1) ?? Q.zero).add(base.get(-1) ?? Q.zero);
    const marg = new Map<string, Q>();
    for (const n of c.model.nodes) if (n.var) marg.set(n.name, (base.get(n.var) ?? Q.zero).div(z));
    const ev = new Map([[c.model.nodes[c.model.nodes.length - 1].name, true]]);
    for (const dot of [modelDot(c.model, marg, marg, new Map()), modelDot(c.model, null, null, new Map()),
      modelDot(c.model, marg, marg, ev, true)])
      assertNodesLabelled(renderSvg(dot), c.model);
  });

test("errors come back with line and column", async () => {
  const t = await translate("event a 0.5\nb = a & c\nevent d 2\n");
  assert.equal(t.r.status, 1);
  assert.deepEqual(t.model.errors.map((e) => [e.line, e.col]), [[2, 9], [3, 9]]);
  assert.match(t.model.errors[0].msg, /unknown name 'c'/);
});

test("ganak-model.js starts fast enough to run on every edit", async () => {
  await translate(HYD);
  const t0 = performance.now();
  for (let i = 0; i < 3; i++) await translate(HYD);
  const ms = (performance.now() - t0) / 3;
  assert.ok(ms < 1500, `${ms.toFixed(0)} ms per call`);
  console.log(`# ganak-model.js: ${ms.toFixed(0)} ms per call`);
});

test("highlighting: classes, lossless, no errors on the examples", () => {
  const src = 'event e "E \\" x" 2e-5 # c\nd = !e & (x.v | true) -> f; g = atleast 2 of (a, b)\nq = 1/3 @ "open';
  const toks = gmTokens(src).map((t) => [t.cls, src.slice(t.b, t.e)]);
  assert.deepEqual(toks, [
    ["kw", "event"], ["def", "e"], ["str", '"E \\" x"'], ["num", "2e-5"], ["com", "# c"],
    ["def", "d"], ["op", "="], ["op", "!"], ["name", "e"], ["op", "&"], ["punc", "("], ["name", "x.v"], ["op", "|"],
    ["const", "true"], ["punc", ")"], ["op", "->"], ["name", "f"], ["punc", ";"],
    ["def", "g"], ["op", "="], ["kw2", "atleast"], ["num", "2"], ["kw2", "of"], ["punc", "("], ["name", "a"],
    ["punc", ","], ["name", "b"], ["punc", ")"],
    ["def", "q"], ["op", "="], ["num", "1/3"], ["bad", "@"], ["bad", '"open'],
  ]);
  const unesc = (h: string) => h.replace(/<[^>]*>/g, "").replace(/&lt;/g, "<").replace(/&gt;/g, ">").replace(/&amp;/g, "&");
  for (const x of MODEL_EXAMPLES) {
    assert.equal(unesc(gmHighlight(x.text)), x.text + "\n", x.name);
    assert.ok(!gmTokens(x.text).some((t) => t.cls === "bad"), x.name);
  }
});
