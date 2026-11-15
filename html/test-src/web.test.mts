// Tests of the web page's logic and of the wasm tools it runs, in Node: the
// page's examples, parsers and query definitions (app.js) against ganak.js and
// bdnnf.js, run the way the page's worker runs them. `npm test` in html/.

import { test } from "node:test";
import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import { createRequire } from "node:module";
import path from "node:path";
import { fileURLToPath } from "node:url";
import { EXAMPLES, QUERIES, parseCnf, fracDecimal, parseModels } from "../app.js";

const HTML = path.dirname(path.dirname(fileURLToPath(import.meta.url)));
const req = createRequire(import.meta.url);

interface Run {
  status: number;
  stdout: string;
  stderr: string;
  files: Record<string, Uint8Array>;
}

// One Emscripten CLI module per call, as the page's worker does.
function runTool(tool: string, args: string[], files: Record<string, string | Uint8Array>,
    outFiles: string[] = []): Promise<Run> {
  return new Promise((resolve) => {
    const js = path.join(HTML, tool + ".js");
    const out: string[] = [];
    const err: string[] = [];
    const mod: any = {
      noInitialRun: true,
      arguments: [],
      print: (t: string) => out.push(t),
      printErr: (t: string) => err.push(t),
      quit: (_s: number, e: unknown) => { throw e; },
      onRuntimeInitialized: () => {
        for (const [p, d] of Object.entries(files)) mod.FS.writeFile(p, d);
        let status = 0;
        try {
          const r = mod.callMain([...args]);  // callMain changes the array
          if (typeof r === "number") status = r;
        } catch (x: any) {
          if (x && typeof x.status === "number") status = x.status;
          else { err.push(String(x)); status = -1; }
        }
        process.exitCode = 0;  // Emscripten sets it on exit(); the tests decide
        const got: Record<string, Uint8Array> = {};
        for (const p of outFiles) { try { got[p] = mod.FS.readFile(p); } catch { /* not written */ } }
        resolve({ status, stdout: out.join("\n"), stderr: err.join("\n"), files: got });
      },
    };
    const argv = process.argv;
    process.argv = [argv[0]];
    new Function("Module", "require", "__dirname", "__filename", readFileSync(js, "utf8"))(mod, req, HTML, js);
    process.argv = argv;
  });
}

const example = (prefix: string) => {
  const e = EXAMPLES.find((x) => x.name.startsWith(prefix));
  assert.ok(e, prefix);
  return e!.cnf;
};
const CAR = example("Car");
const PWMC = example("Weighted");
const circuits = new Map<string, Uint8Array>();

async function compile(cnf: string) {
  const r = await runTool("ganak", ["--verb", "0", "--mode", "1", "--compile", "/out.bdnnf", "/in.cnf"],
    { "/in.cnf": cnf }, ["/out.bdnnf"]);
  assert.equal(r.status, 0, r.stderr);
  const m = r.stdout.match(/c s exact arb \S+ (\S+)/);
  return { count: m ? m[1] : "", bdnnf: r.files["/out.bdnnf"] };
}

async function circuit(cnf: string) {
  if (!circuits.has(cnf)) circuits.set(cnf, (await compile(cnf)).bdnnf);
  return circuits.get(cnf)!;
}

async function bdnnf(cnf: string, args: string[], extra: Record<string, string> = {}, outFiles: string[] = []) {
  return runTool("bdnnf", args, { "/c.bdnnf": await circuit(cnf), "/in.cnf": cnf, ...extra }, outFiles);
}

// the page's own command for query `id`, with its default parameters
function pageArgs(id: string) {
  const q = QUERIES.find((x) => x.id === id)!;
  const p: Record<string, string> = {};
  for (const x of q.params) p[x.key] = x.value;
  return { q, p, args: q.args(p), files: q.files ? q.files(p) : {} };
}

function clauses(cnf: string) {
  return cnf.split("\n").filter((l) => l.trim() && !/^[cp]/.test(l))
    .map((l) => l.trim().split(/\s+/).map(Number).filter((x) => x !== 0));
}
const satisfies = (m: number[], cls: number[][]) => cls.every((c) => c.some((l) => m.includes(l)));

test("parseCnf reads variables, names, weights and projection", () => {
  const car = parseCnf(CAR);
  assert.equal(car.nvars, 12);
  assert.equal(car.names.get(3), "Electric");
  assert.equal(car.weighted, false);
  assert.equal(car.projection, null);
  const w = parseCnf(PWMC);
  assert.equal(w.weighted, true);
  assert.deepEqual(w.projection, [1, 2, 3]);
});

test("fracDecimal approximates fractions", () => {
  assert.equal(fracDecimal("1/2"), "≈ 0.50000000");
  assert.equal(fracDecimal("45"), "");
});

test("parseModels reads DIMACS model lines", () => {
  assert.deepEqual(parseModels("1 -2 3 0\nc comment\n-1 2 0\n"), [[1, -2, 3], [-1, 2]]);
});

const counts: [string, string][] = [["Car", "45"], ["Weighted", "62/125"], ["Parity", "128"],
  ["Graph", "30"], ["Pigeonhole", "0"]];
for (const [name, want] of counts)
  test(`example "${name}" compiles, with Ganak's count ${want}`, async () => {
    const r = await compile(example(name));
    assert.equal(r.count, want);
    assert.ok(r.bdnnf && r.bdnnf.length > 0);
  });

test("every page query runs on the car example", async () => {
  for (const q of QUERIES) {
    const { args, files } = pageArgs(q.id);
    const r = await bdnnf(CAR, args, files, q.id === "graph" ? ["/c.dot"] : []);
    assert.equal(r.status, 0, `${q.id}: ${r.stderr}`);
  }
});

test("verify accepts the compiled circuit", async () => {
  const r = await bdnnf(CAR, ["verify", "/c.bdnnf", "-q"]);
  assert.equal(r.status, 0, r.stdout + r.stderr);
});

test("count, and count under an assumption", async () => {
  assert.equal((await bdnnf(CAR, ["count", "/c.bdnnf"])).stdout.trim(), "45");
  assert.equal((await bdnnf(CAR, ["count", "/c.bdnnf", "--assume", "3"])).stdout.trim(), "9");
});

test("weighted, projected count with the CNF's weights", async () => {
  assert.equal((await bdnnf(PWMC, ["count", "/c.bdnnf", "--cnf", "/in.cnf"])).stdout.trim(), "62/125");
});

test("sample: the page's 20 models all satisfy the CNF", async () => {
  const { args } = pageArgs("sample");
  const ms = parseModels((await bdnnf(CAR, args)).stdout);
  assert.equal(ms.length, 20);
  for (const m of ms) assert.ok(satisfies(m, clauses(CAR)), m.join(" "));
});

test("enumerate lists the 45 models, distinct", async () => {
  const ms = parseModels((await bdnnf(CAR, ["enum", "/c.bdnnf", "-n", "100"])).stdout);
  assert.equal(ms.length, 45);
  assert.equal(new Set(ms.map((m) => m.join(" "))).size, 45);
});

test("feature counts: Electric is in 9 models, Car in all 45", async () => {
  const c = new Map((await bdnnf(CAR, ["lits", "/c.bdnnf"])).stdout.trim().split("\n")
    .map((l) => l.split(/\s+/) as [string, string]));
  assert.equal(c.get("3"), "9");
  assert.equal(c.get("1"), "45");
  assert.equal(c.get("-1"), "0");
});

test("core: Car, Engine and Transmission; no dead feature", async () => {
  const ls = (await bdnnf(CAR, ["core", "/c.bdnnf"])).stdout.split(/\s+/).map(Number).filter((x) => x !== 0);
  assert.deepEqual(ls.sort((a, b) => a - b), [1, 2, 6]);
});

test("atomic sets: the core features form one", async () => {
  const sets = (await bdnnf(CAR, ["atomic", "/c.bdnnf"])).stdout.trim().split("\n")
    .map((l) => l.split(/\s+/).map(Number).filter((x) => x !== 0));
  assert.ok(sets.some((s) => [1, 2, 6].every((v) => s.includes(v))), JSON.stringify(sets));
});

test("query batch: the page's queries, counts and SAT", async () => {
  const { args, files } = pageArgs("batch");
  assert.deepEqual((await bdnnf(CAR, args, files)).stdout.trim().split(/\s+/), ["0", "9", "0"]);
  const sat = await bdnnf(CAR, [...args, "--sat"], files);
  assert.deepEqual(sat.stdout.trim().split(/\s+/), ["UNSAT", "SAT", "UNSAT"], sat.stderr);
});

test("optimise: k best models, values descending and equal to the sum", async () => {
  const { args, files } = pageArgs("opt");
  const obj = new Map(files["/obj.txt"].trim().split("\n").map((l) => l.split(/\s+/).map(Number) as [number, number]));
  const ls = (await bdnnf(CAR, args, files)).stdout.trim().split("\n");
  const vals: number[] = [];
  for (let i = 0; i + 1 < ls.length; i += 2) {
    const v = Number(ls[i].replace("c value ", ""));
    const m = ls[i + 1].split(/\s+/).map(Number).filter((x) => x !== 0);
    assert.equal(m.reduce((s, l) => s + (obj.get(l) ?? 0), 0), v);
    vals.push(v);
  }
  assert.equal(vals.length, 3);
  assert.deepEqual([...vals].sort((a, b) => b - a), vals);
});

test("t-wise: samples cover every satisfiable pair (twise-check)", async () => {
  const r = await bdnnf(CAR, ["twise", "/c.bdnnf", "-t", "2", "-o", "/tw.txt"], {}, ["/tw.txt"]);
  assert.equal(r.status, 0, r.stderr);
  const chk = await bdnnf(CAR, ["twise-check", "/c.bdnnf", "/tw.txt", "-t", "2"],
    { "/tw.txt": new TextDecoder().decode(r.files["/tw.txt"]) });
  assert.equal(chk.status, 0, chk.stdout);
});

test("draw circuit writes a Graphviz graph", async () => {
  const r = await bdnnf(CAR, ["dot", "/c.bdnnf", "/c.dot"], {}, ["/c.dot"]);
  assert.match(new TextDecoder().decode(r.files["/c.dot"]), /^digraph/m);
});

test("circuit text names the format and the variables", async () => {
  assert.match((await bdnnf(CAR, ["dump", "/c.bdnnf"])).stdout, /^bdnnf 1 nvars 12/);
});

test("an unsatisfiable formula: count 0, sampling explains why it cannot", async () => {
  const ph = example("Pigeonhole");
  assert.equal((await bdnnf(ph, ["count", "/c.bdnnf"])).stdout.trim(), "0");
  const r = await bdnnf(ph, ["sample", "/c.bdnnf", "-n", "5"]);
  assert.notEqual(r.status, 0);
  assert.match(r.stderr, /no model/);
});

test("without the cache (--cache 0) the answers are the same", async () => {
  const a = await bdnnf(CAR, ["lits", "/c.bdnnf"]);
  const b = await bdnnf(CAR, ["lits", "/c.bdnnf", "--cache", "0"]);
  assert.equal(a.stdout, b.stdout);
});
