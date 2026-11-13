// Ganak in the browser: count a CNF, compile it to a BDNNF d-DNNF circuit, and
// query the circuit with bdnnf. Every tool run is one fresh worker (worker.ts).

interface RunResult {
  status: number;
  stdout: string;
  stderr: string;
  files: Record<string, Uint8Array>;
  ms: number;
}

const EXAMPLES: { name: string; cnf: string }[] = [
  {
    name: "Car configurator (feature model)",
    cnf: `c A small feature model: a car with an engine, a transmission and options.
c 1 Car
c 2 Engine
c 3 Electric
c 4 Petrol
c 5 Hybrid
c 6 Transmission
c 7 Manual
c 8 Automatic
c 9 Navigation
c 10 Sunroof
c 11 Cabrio
c 12 Screen
p cnf 12 22
1 0
-1 2 0
-2 1 0
-3 2 0
-4 2 0
-5 2 0
-2 3 4 5 0
-3 -4 0
-3 -5 0
-4 -5 0
-1 6 0
-6 1 0
-7 6 0
-8 6 0
-6 7 8 0
-7 -8 0
-9 1 0
-10 1 0
-11 1 0
-12 1 0
-3 8 0
-10 -11 0
-9 12 0
`,
  },
  {
    name: "Weighted, projected (pwmc)",
    cnf: `c t pwmc
p cnf 6 5
1 2 0
-3 -1 -2 0
1 4 5 0
-4 5 0
1 3 0
c p show 1 2 3 0
c p weight 1 4/10 0
c p weight 2 6/10 0
c p weight 3 8/10 0
`,
  },
  {
    name: "Parity chain (XOR of 8 vars is 1)",
    cnf: `c x1 xor ... xor x8 = 1, via t_i = t_{i-1} xor x_{i+1} (vars 9..15)
p cnf 15 29
1 2 -9 0
1 -2 9 0
-1 2 9 0
-1 -2 -9 0
9 3 -10 0
9 -3 10 0
-9 3 10 0
-9 -3 -10 0
10 4 -11 0
10 -4 11 0
-10 4 11 0
-10 -4 -11 0
11 5 -12 0
11 -5 12 0
-11 5 12 0
-11 -5 -12 0
12 6 -13 0
12 -6 13 0
-12 6 13 0
-12 -6 -13 0
13 7 -14 0
13 -7 14 0
-13 7 14 0
-13 -7 -14 0
14 8 -15 0
14 -8 15 0
-14 8 15 0
-14 -8 -15 0
15 0
c p show 1 2 3 4 5 6 7 8 0
`,
  },
  {
    name: "Graph 3-colouring (a 5-cycle)",
    cnf: `c vertex v in {1..5} gets colour c in {1,2,3}: var 3(v-1)+c
p cnf 15 35
1 2 3 0
4 5 6 0
7 8 9 0
10 11 12 0
13 14 15 0
-1 -2 0
-1 -3 0
-2 -3 0
-4 -5 0
-4 -6 0
-5 -6 0
-7 -8 0
-7 -9 0
-8 -9 0
-10 -11 0
-10 -12 0
-11 -12 0
-13 -14 0
-13 -15 0
-14 -15 0
-1 -4 0
-2 -5 0
-3 -6 0
-4 -7 0
-5 -8 0
-6 -9 0
-7 -10 0
-8 -11 0
-9 -12 0
-10 -13 0
-11 -14 0
-12 -15 0
-13 -1 0
-14 -2 0
-15 -3 0
`,
  },
  {
    name: "Pigeonhole: 4 pigeons, 3 holes (UNSAT)",
    cnf: `c pigeon p in hole h: var 3(p-1)+h
p cnf 12 22
1 2 3 0
4 5 6 0
7 8 9 0
10 11 12 0
-1 -4 0
-1 -7 0
-1 -10 0
-4 -7 0
-4 -10 0
-7 -10 0
-2 -5 0
-2 -8 0
-2 -11 0
-5 -8 0
-5 -11 0
-8 -11 0
-3 -6 0
-3 -9 0
-3 -12 0
-6 -9 0
-6 -12 0
-9 -12 0
`,
  },
];

const $ = <T extends HTMLElement>(id: string) => document.getElementById(id) as T;

function el<K extends keyof HTMLElementTagNameMap>(tag: K, attrs: Record<string, string> = {},
    ...kids: (Node | string)[]): HTMLElementTagNameMap[K] {
  const e = document.createElement(tag);
  for (const [k, v] of Object.entries(attrs)) {
    if (k === "class") e.className = v;
    else e.setAttribute(k, v);
  }
  for (const k of kids) e.append(k);
  return e;
}

function runTool(tool: string, args: string[], files: Record<string, string | Uint8Array>,
    outFiles: string[] = []): Promise<RunResult> {
  return new Promise((resolve, reject) => {
    const w = new Worker("worker.js");
    const t0 = performance.now();
    w.onmessage = (e) => {
      w.terminate();
      resolve({ ...e.data, ms: performance.now() - t0 });
    };
    w.onerror = (e) => {
      w.terminate();
      reject(new Error(e.message || "worker error"));
    };
    w.postMessage({ tool, args, files, outFiles });
  });
}

// ---- the CNF ----

interface CnfInfo {
  nvars: number;
  names: Map<number, string>;
  weighted: boolean;
  projection: number[] | null;
}

function parseCnf(text: string): CnfInfo {
  let nvars = 0;
  const names = new Map<number, string>();
  let weighted = false;
  let projection: number[] | null = null;
  for (const raw of text.split("\n")) {
    const t = raw.trim().split(/\s+/);
    if (t[0] === "p" && t.length >= 3) nvars = parseInt(t[2]);
    else if (t[0] === "c" && t[1] === "p" && t[2] === "weight") weighted = true;
    else if (t[0] === "c" && t[1] === "p" && t[2] === "show")
      projection = (projection ?? ([] as number[])).concat(t.slice(3).map(Number).filter((x) => x > 0));
    else if (t[0] === "c" && t.length === 3 && /^\d+$/.test(t[1]) && !/^\d+$/.test(t[2]))
      names.set(parseInt(t[1]), t[2]);
  }
  return { nvars, names, weighted, projection };
}

// ---- state ----

let cnfText = "";
let info: CnfInfo = parseCnf("");
let circuit: Uint8Array | null = null;
let busy = false;

function litName(l: number): string {
  const n = info.names.get(Math.abs(l));
  const v = n ?? String(Math.abs(l));
  return l < 0 ? "¬" + v : v;
}

function setBusy(b: boolean, msg = "") {
  busy = b;
  document.querySelectorAll<HTMLButtonElement>("button.needs-idle").forEach((x) => (x.disabled = b));
  document.querySelectorAll<HTMLButtonElement>("button.needs-circuit").forEach((x) => (x.disabled = b || !circuit));
  $("status").textContent = msg;
}

function log(title: string, r: RunResult) {
  const pre = el("pre", {}, `$ ${title}  (${r.ms.toFixed(0)} ms, exit ${r.status})\n` +
    (r.stdout ? r.stdout + "\n" : "") + (r.stderr ? r.stderr + "\n" : ""));
  $("log").prepend(pre);
}

function showResult(title: string, ms: number, body: Node) {
  const box = $("result");
  box.innerHTML = "";
  box.append(el("div", { class: "result-head" }, el("strong", {}, title), el("span", { class: "muted" }, ` ${ms.toFixed(0)} ms`)), body);
}

function showError(title: string, r: RunResult | null, msg = "") {
  const lines = (r?.stderr ?? "").split("\n").filter((l) => l && !l.startsWith("c o "));
  showResult(title, r?.ms ?? 0, el("div", { class: "error" }, msg || lines.join("\n") || "failed"));
}

// ---- rendering ----

function fracDecimal(s: string): string {
  const [n, d] = s.split("/");
  if (!d) return "";
  const N = BigInt(n), D = BigInt(d);
  const scale = 10n ** 12n;
  const q = (N * scale) / D;
  const x = Number(q) / 1e12;
  return "≈ " + (Math.abs(x) < 1e-6 && N !== 0n ? (Number(N) / Number(D)).toExponential(6) : x.toPrecision(8));
}

function bigNumber(s: string): HTMLElement {
  const box = el("div", { class: "big" }, s.length > 60 ? s.slice(0, 30) + "…" + s.slice(-20) : s);
  const dec = fracDecimal(s);
  return el("div", {}, box, dec ? el("div", { class: "muted" }, dec) : "",
    s.length > 60 ? el("div", { class: "muted" }, `${s.replace("/", "").length} digits`) : "");
}

function parseModels(out: string): number[][] {
  return out.split("\n").filter((l) => l.trim() && !l.startsWith("c")).map((l) =>
    l.trim().split(/\s+/).map(Number).filter((x) => x !== 0));
}

function modelGrid(models: number[][], extra: string[] = []): HTMLElement {
  if (!models.length) return el("div", { class: "muted" }, "no models");
  const vars = Array.from(new Set(models.flatMap((m) => m.map(Math.abs)))).sort((a, b) => a - b);
  const head = el("tr", {}, el("th", {}, "#"), ...(extra.length ? [el("th", {}, "value")] : []),
    ...vars.map((v) => el("th", { title: String(v) }, info.names.get(v) ?? String(v))));
  const shown = models.slice(0, 500);
  const rows = shown.map((m, i) => {
    const set = new Set(m);
    return el("tr", {}, el("td", { class: "muted" }, String(i + 1)), ...(extra.length ? [el("td", {}, extra[i])] : []),
      ...vars.map((v) => el("td", { class: set.has(v) ? "on" : "off", title: litName(set.has(v) ? v : -v) },
        set.has(v) ? "1" : "0")));
  });
  const table = el("table", { class: "grid" }, el("thead", {}, head), el("tbody", {}, ...rows));
  return el("div", { class: "scroll" }, table,
    models.length > shown.length ? el("div", { class: "muted" }, `first ${shown.length} of ${models.length} shown`) : "");
}

function chips(lits: number[], cls = ""): HTMLElement {
  if (!lits.length) return el("div", { class: "muted" }, "none");
  return el("div", { class: "chips" }, ...lits.map((l) => el("span", { class: "chip " + cls + (l < 0 ? " neg" : "") }, litName(l))));
}

// ---- compile and count ----

async function countWithGanak() {
  setBusy(true, "Counting...");
  try {
    const r = await runTool("ganak", ["--verb", "1", "--mode", "1", "/in.cnf"], { "/in.cnf": cnfText });
    log("ganak --mode 1 in.cnf", r);
    const m = r.stdout.match(/c s exact arb \S+ (\S+)/);
    if (r.status !== 0 || !m) showError("Count", r);
    else showResult("Model count (Ganak)", r.ms, bigNumber(m[1]));
  } finally {
    setBusy(false);
  }
}

async function compile() {
  setBusy(true, "Compiling...");
  circuit = null;
  $("circuit-info").innerHTML = "";
  try {
    const r = await runTool("ganak", ["--verb", "1", "--mode", "1", "--compile", "/out.bdnnf", "/in.cnf"],
      { "/in.cnf": cnfText }, ["/out.bdnnf"]);
    log("ganak --mode 1 --compile out.bdnnf in.cnf", r);
    const m = r.stdout.match(/c s exact arb \S+ (\S+)/);
    if (r.status !== 0 || !r.files["/out.bdnnf"]) {
      showError("Compile", r);
      return;
    }
    circuit = r.files["/out.bdnnf"];
    $("circuit-hint").textContent = "";
    const v = await bdnnf(["verify", "/c.bdnnf", "-q"], false);
    const d = await bdnnf(["dump", "/c.bdnnf"], false);
    const short: Record<string, string> = { DECOMPOSABLE: "DEC", DETERMINISTIC: "DET", SMOOTH: "SMOOTH", CLEAN: "CLEAN" };
    const props = (d?.stdout.match(/^props (.*)$/m) ?? ["", ""])[1].split(/\s+/).map((p) => short[p] ?? p).join(" ");
    const nodes = d ? d.stdout.split("\n").filter((l) => /^\d+ (and|or)/.test(l)).length : 0;
    $("circuit-info").append(
      el("div", { class: "stats" },
        stat("count", m ? m[1] : "?"), stat("compile", `${r.ms.toFixed(0)} ms`),
        stat("nodes", String(nodes)), stat("size", `${circuit.length} B`),
        stat("properties", props || "-"), stat("verify", v && v.status === 0 ? "OK" : "FAILED")));
    showResult("Compiled", r.ms, el("div", {}, "The circuit is ready: pick a query on the right."));
  } finally {
    setBusy(false);
  }
}

function stat(k: string, v: string): HTMLElement {
  return el("div", { class: "stat" }, el("div", { class: "muted" }, k),
    el("div", { class: "stat-v", title: v }, v.length > 24 ? v.slice(0, 22) + "…" : v));
}

async function bdnnf(args: string[], show = true, extraFiles: Record<string, string> = {},
    outFiles: string[] = []): Promise<RunResult | null> {
  if (!circuit) return null;
  const r = await runTool("bdnnf", args, { "/c.bdnnf": circuit, "/in.cnf": cnfText, ...extraFiles }, outFiles);
  if (show) log("bdnnf " + args.join(" "), r);
  return r;
}

// ---- queries ----

interface Param {
  key: string;
  label: string;
  kind: "number" | "text" | "check" | "area";
  value: string;
  hint?: string;
}

interface Query {
  id: string;
  label: string;
  desc: string;
  params: Param[];
  assume: boolean;
  args: (p: Record<string, string>) => string[];
  render: (r: RunResult, p: Record<string, string>) => Node;
  files?: (p: Record<string, string>) => Record<string, string>;
}

function lines(s: string) {
  return s.split("\n").map((l) => l.trim()).filter((l) => l);
}

const QUERIES: Query[] = [
  {
    id: "count", label: "Count", desc: "The (weighted, projected) model count.", params: [], assume: true,
    args: () => ["count", "/c.bdnnf"],
    render: (r) => bigNumber(lines(r.stdout)[0] ?? "?"),
  },
  {
    id: "sample", label: "Sample", desc: "Random models, uniformly (or by weight when weights are used).",
    params: [{ key: "n", label: "models", kind: "number", value: "20" }, { key: "seed", label: "seed", kind: "number", value: "1" }],
    assume: true, args: (p) => ["sample", "/c.bdnnf", "-n", p.n, "--seed", p.seed],
    render: (r) => modelGrid(parseModels(r.stdout)),
  },
  {
    id: "enum", label: "Enumerate", desc: "The models in a fixed order, from model number “start”.",
    params: [{ key: "n", label: "models", kind: "number", value: "50" }, { key: "start", label: "start", kind: "number", value: "0" }],
    assume: true, args: (p) => ["enum", "/c.bdnnf", "-n", p.n, "--start", p.start],
    render: (r) => modelGrid(parseModels(r.stdout)),
  },
  {
    id: "lits", label: "Feature counts", desc: "For each literal, the count of the models containing it.", params: [], assume: true,
    args: () => ["lits", "/c.bdnnf"],
    render: (r) => {
      const c = new Map<number, string>();
      for (const l of lines(r.stdout)) {
        const [a, b] = l.split(/\s+/);
        c.set(Number(a), b);
      }
      const vars = Array.from(c.keys()).filter((x) => x > 0).sort((a, b) => a - b);
      const toNum = (s: string) => {
        const [n, d] = s.split("/");
        return d ? Number(BigInt(n) * 1000000n / BigInt(d)) / 1e6 : Number(n);
      };
      const rows = vars.map((v) => {
        const p = toNum(c.get(v)!), q = toNum(c.get(-v)!);
        const share = p + q > 0 ? p / (p + q) : 0;
        return el("tr", {}, el("td", {}, litName(v)), el("td", { class: "num" }, c.get(v)!), el("td", { class: "num" }, c.get(-v)!),
          el("td", {}, el("div", { class: "bar" }, el("div", { style: `width:${(share * 100).toFixed(1)}%` })),
            el("span", { class: "muted" }, ` ${(share * 100).toFixed(1)}%`)));
      });
      return el("div", { class: "scroll" }, el("table", { class: "list" },
        el("thead", {}, el("tr", {}, el("th", {}, "feature"), el("th", {}, "with"), el("th", {}, "without"), el("th", {}, "share"))),
        el("tbody", {}, ...rows)));
    },
  },
  {
    id: "core", label: "Core / dead", desc: "Features in every model (core) and in none (dead).", params: [], assume: true,
    args: () => ["core", "/c.bdnnf"],
    render: (r) => {
      const ls = (lines(r.stdout)[0] ?? "").split(/\s+/).map(Number).filter((x) => x !== 0);
      return el("div", {}, el("h4", {}, "core"), chips(ls.filter((l) => l > 0), "core"),
        el("h4", {}, "dead"), chips(ls.filter((l) => l < 0).map((l) => -l), "dead"));
    },
  },
  {
    id: "atomic", label: "Atomic sets", desc: "Groups of features that are equal in every model.",
    params: [{ key: "cross", label: "also opposite (--cross)", kind: "check", value: "" }], assume: true,
    args: (p) => ["atomic", "/c.bdnnf", ...(p.cross ? ["--cross"] : [])],
    render: (r) => {
      const sets = lines(r.stdout).map((l) => l.split(/\s+/).map(Number).filter((x) => x !== 0));
      return sets.length ? el("div", {}, ...sets.map((s) => chips(s, "atom"))) : el("div", { class: "muted" }, "no atomic sets");
    },
  },
  {
    id: "batch", label: "Query batch", desc: "One count (or SAT/UNSAT) per line of literals.",
    params: [{ key: "q", label: "queries", kind: "area", value: "1 -2\n3\n-1 -3" }, { key: "sat", label: "SAT only", kind: "check", value: "" }],
    assume: false, files: (p) => ({ "/q.txt": p.q + "\n" }),
    args: (p) => ["batch", "/c.bdnnf", "/q.txt", ...(p.sat ? ["--sat"] : [])],
    render: (r, p) => {
      const qs = lines(p.q);
      const ans = lines(r.stdout);
      return el("table", { class: "list" }, el("thead", {}, el("tr", {}, el("th", {}, "query"), el("th", {}, "answer"))),
        el("tbody", {}, ...qs.map((q, i) => el("tr", {},
          el("td", {}, q.split(/\s+/).map(Number).map(litName).join(" ∧ ")), el("td", { class: "num" }, ans[i] ?? "")))));
    },
  },
  {
    id: "opt", label: "Optimise", desc: "The best models under a sum of per-literal values (lines of “literal value”).",
    params: [{ key: "obj", label: "objective", kind: "area", value: "3 5\n4 2\n9 3\n-11 1" },
      { key: "k", label: "best k", kind: "number", value: "3" }, { key: "min", label: "minimise", kind: "check", value: "" }],
    assume: true, files: (p) => ({ "/obj.txt": p.obj + "\n" }),
    args: (p) => ["opt", "/c.bdnnf", "--obj", "/obj.txt", "-k", p.k, ...(p.min ? ["--min"] : [])],
    render: (r) => {
      const ls = lines(r.stdout);
      const vals: string[] = [];
      const ms: number[][] = [];
      for (let i = 0; i + 1 < ls.length; i += 2) {
        vals.push(ls[i].replace(/^c value /, ""));
        ms.push(ls[i + 1].split(/\s+/).map(Number).filter((x) => x !== 0));
      }
      return modelGrid(ms, vals);
    },
  },
  {
    id: "twise", label: "t-wise sample", desc: "Samples covering every satisfiable t-tuple of literals.",
    params: [{ key: "t", label: "t", kind: "number", value: "2" }, { key: "seed", label: "seed", kind: "number", value: "1" }],
    assume: true, args: (p) => ["twise", "/c.bdnnf", "-t", p.t, "--seed", p.seed],
    render: (r) => {
      const ms = parseModels(r.stdout);
      const note = (r.stderr.split("\n").find((l) => l.startsWith("bdnnf twise")) ?? "").replace("bdnnf twise: ", "");
      return el("div", {}, el("div", { class: "muted" }, note), modelGrid(ms));
    },
  },
  {
    id: "graph", label: "Draw circuit", desc: "The circuit as a graph (Graphviz).", params: [], assume: false,
    args: () => ["dot", "/c.bdnnf", "/c.dot"],
    render: () => el("div", {}),
  },
  {
    id: "dump", label: "Circuit text", desc: "The circuit in BDNNF's text form.", params: [], assume: false,
    args: () => ["dump", "/c.bdnnf"],
    render: (r) => el("pre", { class: "dump" }, r.stdout),
  },
];

let current: Query = QUERIES[0];

function buildQueryButtons() {
  const box = $("query-buttons");
  for (const q of QUERIES) {
    const b = el("button", { class: "qbtn needs-circuit", title: q.desc }, q.label);
    b.addEventListener("click", () => selectQuery(q));
    box.append(b);
  }
}

function selectQuery(q: Query) {
  current = q;
  document.querySelectorAll(".qbtn").forEach((b) => b.classList.toggle("active", b.textContent === q.label));
  const form = $("query-form");
  form.innerHTML = "";
  form.append(el("p", { class: "muted" }, q.desc));
  for (const p of q.params) {
    const id = `p-${q.id}-${p.key}`;
    let input: HTMLElement;
    if (p.kind === "area") input = el("textarea", { id, rows: "4" }, p.value);
    else if (p.kind === "check") input = el("input", { id, type: "checkbox", ...(p.value ? { checked: "" } : {}) });
    else input = el("input", { id, type: p.kind === "number" ? "number" : "text", value: p.value });
    form.append(el("label", { for: id, class: p.kind === "check" ? "inline" : "" }, p.label), input);
  }
  const run = el("button", { class: "primary needs-circuit" }, "Run");
  run.addEventListener("click", runQuery);
  form.append(run);
  run.disabled = busy || !circuit;
}

function readParams(q: Query): Record<string, string> {
  const p: Record<string, string> = {};
  for (const x of q.params) {
    const e = document.getElementById(`p-${q.id}-${x.key}`) as HTMLInputElement | HTMLTextAreaElement;
    p[x.key] = x.kind === "check" ? ((e as HTMLInputElement).checked ? "1" : "") : e.value;
  }
  return p;
}

async function runQuery() {
  const q = current;
  const p = readParams(q);
  const args = q.args(p);
  if (q.assume) {
    const a = $<HTMLInputElement>("assume").value.trim().split(/[\s,]+/).filter((x) => x);
    if (a.length) args.push("--assume", a.join(","));
  }
  if ($<HTMLInputElement>("use-cnf").checked && !["graph", "dump"].includes(q.id)) args.push("--cnf", "/in.cnf");
  setBusy(true, `Running ${q.label.toLowerCase()}...`);
  try {
    const outFiles = q.id === "graph" ? ["/c.dot"] : [];
    const r = await bdnnf(args, true, q.files ? q.files(p) : {}, outFiles);
    if (!r) return;
    if (r.status !== 0) {
      showError(q.label, r);
      return;
    }
    if (q.id === "graph") {
      const dot = new TextDecoder().decode(r.files["/c.dot"]);
      showResult(q.label, r.ms, await renderDot(dot));
      return;
    }
    showResult(q.label, r.ms, q.render(r, p));
  } catch (e) {
    showError(q.label, null, String(e));
  } finally {
    setBusy(false);
  }
}

let viz: any = null;
async function renderDot(dot: string): Promise<Node> {
  if (!viz) {
    const m: any = await import("https://cdn.jsdelivr.net/npm/@viz-js/viz@3.7.0/+esm" as any);
    viz = await m.instance();
  }
  return el("div", { class: "graph" }, viz.renderSVGElement(dot));
}

// ---- wiring ----

function setCnf(text: string) {
  cnfText = text;
  info = parseCnf(text);
  circuit = null;
  $("circuit-info").innerHTML = "";
  $("circuit-hint").textContent = "Compile the formula first.";
  const notes: string[] = [`${info.nvars} vars`];
  if (info.projection) notes.push(`projected on ${info.projection.length}`);
  if (info.weighted) notes.push("weighted");
  if (info.names.size) notes.push(`${info.names.size} named features`);
  $("cnf-info").textContent = notes.join(" · ");
  const useCnf = $<HTMLInputElement>("use-cnf");
  useCnf.checked = info.weighted;
  useCnf.disabled = !info.weighted;
  setBusy(false);
}

function init() {
  const sel = $<HTMLSelectElement>("examples");
  EXAMPLES.forEach((x, i) => sel.append(el("option", { value: String(i) }, x.name)));
  const ta = $<HTMLTextAreaElement>("cnf");
  sel.addEventListener("change", () => {
    ta.value = EXAMPLES[Number(sel.value)].cnf;
    setCnf(ta.value);
  });
  ta.addEventListener("input", () => setCnf(ta.value));
  $<HTMLInputElement>("file").addEventListener("change", async (e) => {
    const f = (e.target as HTMLInputElement).files?.[0];
    if (!f) return;
    ta.value = await f.text();
    setCnf(ta.value);
  });
  $("btn-count").addEventListener("click", countWithGanak);
  $("btn-compile").addEventListener("click", compile);
  $("btn-download").addEventListener("click", () => {
    if (!circuit) return;
    const a = el("a", { href: URL.createObjectURL(new Blob([circuit as BlobPart])), download: "circuit.bdnnf" });
    a.click();
  });
  buildQueryButtons();
  ta.value = EXAMPLES[0].cnf;
  setCnf(ta.value);
  selectQuery(QUERIES[0]);
}

init();
