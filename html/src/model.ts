// The "Probability models" tab: a ganak-model (.gm) model is translated to a
// weighted CNF, compiled once by Ganak, and every answer comes from bdnnf on
// that circuit. Pure helpers are exported for the Node tests.

export interface RunResult {
  status: number;
  stdout: string;
  stderr: string;
  files: Record<string, Uint8Array>;
  ms: number;
}

export interface Deps {
  runTool(tool: string, args: string[], files: Record<string, string | Uint8Array>, outFiles?: string[]): Promise<RunResult>;
  renderDot(dot: string): Promise<HTMLElement>;
  log(title: string, r: RunResult): void;
  openInCnf(cnf: string): void;
}

export interface GmNode {
  name: string;
  label: string;
  kind: "event" | "value" | "flip" | "def";
  var: number;
  line: number;
  anon: boolean;
  p?: string;
  choice?: string;
  text?: string;
  deps: { name: string; pol: "+" | "-" | "*" }[];
}

export interface GmStmt {
  kind: "constraint" | "query" | "require";
  text: string;
  line: number;
  lit: number;
  given: number[];
  op?: string;
  bound?: string;
  deps?: { name: string; pol: "+" | "-" | "*" }[];
}

export interface GmModel {
  errors: { line: number; col: number; msg: string }[];
  nvars: number;
  nodes: GmNode[];
  choices: { name: string; label: string; values: string[] }[];
  stmts: GmStmt[];
}

// ---- exact rationals ----

function gcd(a: bigint, b: bigint): bigint {
  if (a < 0n) a = -a;
  if (b < 0n) b = -b;
  while (b) [a, b] = [b, a % b];
  return a;
}

export class Q {
  readonly n: bigint;
  readonly d: bigint;
  constructor(n: bigint, d: bigint = 1n) {
    if (d === 0n) throw new Error("division by zero");
    if (d < 0n) {
      n = -n;
      d = -d;
    }
    const g = gcd(n, d) || 1n;
    this.n = n / g;
    this.d = d / g;
  }
  static parse(s: string): Q {
    const t = s.trim();
    const i = t.indexOf("/");
    return i < 0 ? new Q(BigInt(t)) : new Q(BigInt(t.slice(0, i)), BigInt(t.slice(i + 1)));
  }
  static readonly zero = new Q(0n);
  static readonly one = new Q(1n);
  add(o: Q) { return new Q(this.n * o.d + o.n * this.d, this.d * o.d); }
  sub(o: Q) { return new Q(this.n * o.d - o.n * this.d, this.d * o.d); }
  mul(o: Q) { return new Q(this.n * o.n, this.d * o.d); }
  div(o: Q) { return new Q(this.n * o.d, this.d * o.n); }
  cmp(o: Q) {
    const x = this.n * o.d - o.n * this.d;
    return x < 0n ? -1 : x > 0n ? 1 : 0;
  }
  isZero() { return this.n === 0n; }
  toString() { return this.d === 1n ? String(this.n) : `${this.n}/${this.d}`; }
  toNumber(): number {
    if (this.n === 0n) return 0;
    const neg = this.n < 0n;
    const n = neg ? -this.n : this.n;
    const shift = this.d.toString().length - n.toString().length + 18;
    const scaled = shift >= 0 ? (n * 10n ** BigInt(shift)) / this.d : n / (this.d * 10n ** BigInt(-shift));
    const x = Number(scaled) * Math.pow(10, -shift);
    return neg ? -x : x;
  }
}

const SUP: Record<string, string> = { "-": "⁻", "0": "⁰", "1": "¹", "2": "²", "3": "³", "4": "⁴", "5": "⁵", "6": "⁶", "7": "⁷", "8": "⁸", "9": "⁹" };

export function fmtP(q: Q | null): string {
  if (!q) return "undefined";
  if (q.isZero()) return "0";
  if (q.cmp(Q.one) === 0) return "1";
  const x = q.toNumber();
  if (x > 0.999 && x <= 1) return "1 − " + fmtP(Q.one.sub(q));
  if (Math.abs(x) >= 1e-3) return String(Number(x.toPrecision(4)));
  const [m, e] = x.toExponential(2).split("e");
  return `${m}×10${e.replace("+", "").split("").map((c) => SUP[c] ?? c).join("")}`;
}

// ---- bdnnf output ----

export function parseLits(out: string): Map<number, Q> {
  const m = new Map<number, Q>();
  for (const l of out.split("\n")) {
    const t = l.trim().split(/\s+/);
    if (t.length !== 2 || l.startsWith("c")) continue;
    m.set(Number(t[0]), Q.parse(t[1]));
  }
  return m;
}

export function batchFile(m: GmModel): string {
  const lines: string[] = [];
  for (const s of m.stmts) {
    if (s.kind === "constraint") continue;
    lines.push(s.given.join(" "));
    lines.push([...s.given, s.lit].join(" "));
  }
  return lines.join("\n") + "\n";
}

export interface Answer {
  stmt: GmStmt;
  p: Q | null;
  pass?: boolean;
}

export function answers(m: GmModel, batchOut: string): Answer[] {
  const vals = batchOut.split("\n").map((l) => l.trim()).filter((l) => l && !l.startsWith("c")).map(Q.parse);
  const out: Answer[] = [];
  let i = 0;
  for (const s of m.stmts) {
    if (s.kind === "constraint") continue;
    const den = vals[i++], num = vals[i++];
    const p = !den || den.isZero() ? null : num.div(den);
    const a: Answer = { stmt: s, p };
    if (s.kind === "require") {
      if (!p) a.pass = false;
      else {
        const c = p.cmp(Q.parse(s.bound!));
        a.pass = s.op === "<" ? c < 0 : s.op === "<=" ? c <= 0 : s.op === ">" ? c > 0 : c >= 0;
      }
    }
    out.push(a);
  }
  return out;
}

export function pretty(s: GmStmt): string {
  return "P(" + s.text.replace(" given ", " | ") + ")";
}

export const randomNodes = (m: GmModel) => m.nodes.filter((n) => n.kind !== "def");

export function nodeWeight(n: GmNode, value: boolean): Q {
  const p = Q.parse(n.p ?? "1");
  if (n.kind === "value") return value ? p : Q.one;
  return value ? p : Q.one.sub(p);
}

export interface Impact {
  node: GmNode;
  p: Q;
  ifYes: Q | null;
  ifNo: Q | null;
}

// From W(l ∧ E) (base) and W(l ∧ T ∧ E) (withTop) for the literals of the
// declared events: P(x | E), P(T | x, E) and P(T | ¬x, E), most influential first.
export function impacts(m: GmModel, base: Map<number, Q>, withTop: Map<number, Q>): { pT: Q | null; rows: Impact[] } {
  const rows: Impact[] = [];
  let pT: Q | null = null;
  for (const n of randomNodes(m)) {
    const wx = base.get(n.var) ?? Q.zero, wnx = base.get(-n.var) ?? Q.zero;
    const tx = withTop.get(n.var) ?? Q.zero, tnx = withTop.get(-n.var) ?? Q.zero;
    const z = wx.add(wnx);
    if (z.isZero()) continue;
    pT ??= tx.add(tnx).div(z);
    if (n.anon) continue;
    rows.push({ node: n, p: wx.div(z), ifYes: wx.isZero() ? null : tx.div(wx), ifNo: wnx.isZero() ? null : tnx.div(wnx) });
  }
  const swing = (r: Impact) => {
    const d = (r.ifYes ?? pT ?? Q.zero).sub(r.ifNo ?? pT ?? Q.zero);
    return d.n < 0n ? new Q(-d.n, d.d) : d;
  };
  rows.sort((a, b) => swing(b).cmp(swing(a)));
  return { pT, rows };
}

// ---- the graph ----

function esc(s: string) {
  return s.replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;").replace(/"/g, "&quot;");
}

function fill(p: number) {
  if (p === 0) return "#f1f3f5";
  if (p >= 0.5) return "#bcd4fb";
  if (p >= 1e-2) return "#d8e6fd";
  if (p >= 1e-5) return "#edf3fe";
  return "#ffffff";
}

export function modelDot(m: GmModel, marg: Map<string, Q> | null, prior: Map<string, Q> | null,
    evidence: Map<string, boolean>): string {
  const idx = new Map(m.nodes.map((n, i) => [n.name, i]));
  const asked = new Set(m.stmts.filter((s) => s.kind !== "constraint").map((s) => Math.abs(s.lit)));
  const L: string[] = ['digraph model {', '  rankdir=LR; bgcolor="transparent"; nodesep=0.18; ranksep=0.5;',
    '  node [fontname="Helvetica", fontsize=11, style="filled", color="#7b8494", penwidth=1];',
    '  edge [color="#8a94a6", arrowsize=0.6];'];
  const nodeLine = (n: GmNode, i: number) => {
    const p = marg?.get(n.name) ?? null;
    const pr = prior?.get(n.name) ?? null;
    let ptxt = p ? "P = " + fmtP(p) : n.p ? "p = " + fmtP(Q.parse(n.p)) : "";
    if (p && pr && !pr.isZero() && evidence.size && p.cmp(pr) !== 0) {
      const r = p.toNumber() / pr.toNumber();
      if (p.isZero()) ptxt += "  ▼";
      else if (r >= 1.5 || r <= 1 / 1.5) ptxt += r >= 1 ? `  ▲×${fmtFactor(r)}` : `  ▼÷${fmtFactor(1 / r)}`;
    }
    const title = n.anon ? (n.label || n.name) : (n.label || n.name);
    const sub = n.label && !n.anon ? `<BR/><FONT POINT-SIZE="8" COLOR="#6b7280">${esc(n.name)}</FONT>` : "";
    const label = `<${esc(title)}${sub}<BR/><FONT POINT-SIZE="10">${esc(ptxt)}</FONT>>`;
    const ev = evidence.get(n.name);
    const shape = asked.has(n.var) ? 'shape=box, style="filled", peripheries=2'
      : n.kind === "def" ? 'shape=box, style="filled,rounded"' : n.anon ? "shape=ellipse, fontsize=9" : "shape=ellipse";
    const border = ev === undefined ? "" : `, penwidth=3, color="${ev ? "#16a34a" : "#dc2626"}"`;
    const tip = n.name + (n.text ? " = " + n.text : n.p ? `, probability ${n.p}` : "");
    return `  n${i} [id="gm${i}", ${shape}, label=${label}, fillcolor="${fill(p ? p.toNumber() : n.p ? Q.parse(n.p).toNumber() : 1)}"${border}, tooltip="${esc(tip)}"];`;
  };
  const inChoice = new Set<string>();
  m.choices.forEach((c, ci) => {
    L.push(`  subgraph cluster_${ci} { label=<${esc(c.label || c.name)}>; fontname="Helvetica"; fontsize=10; style="rounded,dashed"; color="#9aa3af";`);
    for (const v of c.values) {
      inChoice.add(v);
      L.push("  " + nodeLine(m.nodes[idx.get(v)!], idx.get(v)!));
    }
    L.push("  }");
  });
  m.nodes.forEach((n, i) => {
    if (!inChoice.has(n.name)) L.push(nodeLine(n, i));
  });
  m.stmts.forEach((s, k) => {
    if (s.kind !== "constraint") return;
    L.push(`  c${k} [id="gc${k}", shape=hexagon, style="filled", fillcolor="#fff4d6", color="#b7791f", fontsize=10, `
      + `label=<${esc(s.text)}>, tooltip="constraint ${esc(s.text)}"];`);
    for (const d of s.deps ?? []) {
      const j = idx.get(d.name);
      if (j !== undefined) L.push(`  c${k} -> n${j} [color="#b7791f", style=${d.pol === "+" ? "solid" : d.pol === "-" ? "dashed" : "dotted"}${d.pol === "-" ? ", arrowhead=odot" : ""}];`);
    }
  });
  m.nodes.forEach((n, i) => {
    for (const d of n.deps) {
      const j = idx.get(d.name);
      if (j === undefined) continue;
      const st = d.pol === "-" ? ' [style=dashed, arrowhead=odot]' : d.pol === "*" ? " [style=dotted]" : "";
      L.push(`  n${i} -> n${j}${st};`);
    }
  });
  L.push("}");
  return L.join("\n");
}

function fmtFactor(r: number) {
  return r >= 100 ? r.toExponential(1).replace("e+", "e") : r >= 10 ? r.toFixed(0) : r.toFixed(1);
}

// ---- examples ----

export const MODEL_EXAMPLES: { name: string; text: string }[] = [
  {
    name: "Weather (constraints)",
    text: `# Constraints rule out impossible weather; every answer is conditioned on them.
event sun       "The sun shines"       0.6
event rain      "It rains"             0.3
event clouds    "Overcast"             0.5
event wind      "Windy"                0.3
event sprinkler "The sprinkler is on"  0.4

constraint !(rain & sun)          # it never rains while the sun shines
constraint rain -> clouds         # no rain without clouds
constraint !(sun & clouds)        # no sunshine when overcast
constraint rain -> !sprinkler     # the sprinkler has a rain sensor

wet "The grass is wet" = rain | sprinkler
picnic "Picnic in the park" = sun & !wind & !wet
umbrella "Take an umbrella" = rain | clouds & wind

query rain
query sun
query rain given wet
query sprinkler given wet
query picnic
query umbrella given !sun
`,
  },
  {
    name: "Airliner flight-control hydraulics",
    text: `# Hydraulic power for the flight controls of a twin-engine airliner, per flight.
# Three independent hydraulic systems (green, yellow, blue). The numbers are
# illustrative, not data of any real aircraft.

event eng1   "Engine 1 shuts down"              2e-5
event eng2   "Engine 2 shuts down"              2e-5
event dual   "Both engines lost (birds, ash)"   1e-7
event apu    "APU fails to start"               5e-2
event edp_g  "Green engine-driven pump fails"   1e-4
event edp_y  "Yellow engine-driven pump fails"  1e-4
event epmp_y "Yellow electric pump fails"       2e-4
event epmp_b "Blue electric pump fails"         2e-4
event rat    "Ram air turbine fails"            1e-3
event ptu    "Power transfer unit fails"        1e-3
event leak_g "Green fluid loss"                 5e-5
event leak_y "Yellow fluid loss"                5e-5
event leak_b "Blue fluid loss"                  5e-5

eng1_out "Engine 1 out" = eng1 | dual
eng2_out "Engine 2 out" = eng2 | dual
ac_power "AC electrical power" = !eng1_out | !eng2_out | !apu

# Pressure sources. The PTU lets a pressurised system drive the other one.
g_src "Green source" = !eng1_out & !edp_g
y_src "Yellow source" = (!eng2_out & !edp_y) | (ac_power & !epmp_y)
green  "Green pressure"  = !leak_g & (g_src | (y_src & !leak_y & !ptu))
yellow "Yellow pressure" = !leak_y & (y_src | (g_src & !leak_g & !ptu))
blue   "Blue pressure"   = !leak_b & ((ac_power & !epmp_b) | !rat)

elevator "Elevator powered" = green | blue
ailerons "Ailerons powered" = green | blue
rudder   "Rudder powered"   = atleast 1 of (green, blue, yellow)
lost "Loss of control" = !(elevator & ailerons & rudder)
land_asap "Land at the nearest airport" = atmost 1 of (green, blue, yellow)

query lost
query land_asap
query lost given eng1_out
query lost given dual
query dual given !green, !yellow

# Catastrophic failure conditions must stay below 1e-9 per flight. This one
# fails: "What matters" below the graph shows why. Powering the elevator and the
# ailerons from all three systems brings it down to about 6e-12.
require lost < 1e-9
`,
  },
];

// ---- the tab ----

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

const dec = (b: Uint8Array) => new TextDecoder().decode(b);

interface State {
  text: string;
  model: GmModel | null;
  cnf: string;
  circuit: Uint8Array | null;
  z: Q;
  prior: Map<string, Q> | null;
  marg: Map<string, Q> | null;
  pe: Q | null;
  evidence: Map<string, boolean>;
  gen: number;
  ev: number;
  wm: number;
  baseLits: Map<number, Q> | null;
  baseEv: number;
}

const S: State = { text: "", model: null, cnf: "", circuit: null, z: Q.one, prior: null, marg: null, pe: null,
  evidence: new Map(), gen: 0, ev: 0, wm: 0, baseLits: null, baseEv: -1 };
let D: Deps;
let timer: number | undefined;

function status(msg: string) {
  $("gm-status").textContent = msg;
}

function pq(q: Q | null, cls = "") {
  return el("span", { class: "p " + cls, title: q ? q.toString() : "" }, fmtP(q));
}

async function bdnnf(args: string[], extra: Record<string, string> = {}, show = true): Promise<RunResult> {
  const r = await D.runTool("bdnnf", args, { "/c.bdnnf": S.circuit!, "/m.cnf": S.cnf, ...extra });
  if (show) D.log("bdnnf " + args.join(" "), r);
  return r;
}

function nodeByName(name: string) {
  return S.model!.nodes.find((n) => n.name === name);
}

function evidenceLits(): number[] {
  const out: number[] = [];
  for (const [name, v] of S.evidence) {
    const n = nodeByName(name);
    if (n) out.push(v ? n.var : -n.var);
  }
  return out;
}

function assumeArgs(lits: number[]) {
  return lits.length ? ["--assume", lits.join(",")] : [];
}

function margMap(lits: Map<number, Q>): { marg: Map<string, Q>; total: Q } {
  const m = S.model!;
  let total = Q.zero;
  const first = m.nodes.find((n) => n.var);
  if (first) total = (lits.get(first.var) ?? Q.zero).add(lits.get(-first.var) ?? Q.zero);
  const marg = new Map<string, Q>();
  for (const n of m.nodes) if (n.var && !total.isZero()) marg.set(n.name, (lits.get(n.var) ?? Q.zero).div(total));
  return { marg, total };
}

function showErrors(errs: GmModel["errors"], extra = "") {
  const box = $("gm-errors");
  box.innerHTML = "";
  if (!errs.length && !extra) {
    box.hidden = true;
    return;
  }
  box.hidden = false;
  if (extra) box.append(el("div", {}, extra));
  for (const e of errs) {
    const a = el("a", { href: "#" }, `line ${e.line}:${e.col}`);
    a.addEventListener("click", (ev) => {
      ev.preventDefault();
      selectLine(e.line);
    });
    box.append(el("div", {}, a, " " + e.msg));
  }
}

function selectLine(line: number) {
  const ta = $<HTMLTextAreaElement>("gm-text");
  const ls = ta.value.split("\n");
  let a = 0;
  for (let i = 0; i < line - 1 && i < ls.length; i++) a += ls[i].length + 1;
  ta.focus();
  ta.setSelectionRange(a, a + (ls[line - 1]?.length ?? 0));
  const lh = parseFloat(getComputedStyle(ta).lineHeight) || 18;
  ta.scrollTop = Math.max(0, (line - 4) * lh);
}

export async function compute() {
  const g = ++S.gen;
  const text = $<HTMLTextAreaElement>("gm-text").value;
  S.text = text;
  status("Translating...");
  const r1 = await D.runTool("ganak-model", ["/m.gm", "-o", "/m.cnf", "--json", "/m.json"], { "/m.gm": text },
    ["/m.cnf", "/m.json"]);
  if (g !== S.gen) return;
  D.log("ganak-model m.gm -o m.cnf --json m.json", r1);
  if (!r1.files["/m.json"]) {
    showErrors([], r1.stderr || "ganak-model failed");
    status("");
    return;
  }
  const model: GmModel = JSON.parse(dec(r1.files["/m.json"]));
  showErrors(model.errors);
  if (model.errors.length || r1.status !== 0 || !r1.files["/m.cnf"]) {
    status(model.errors.length ? `${model.errors.length} error${model.errors.length > 1 ? "s" : ""}` : "");
    $("gm-results").classList.add("stale");
    return;
  }
  const cnf = dec(r1.files["/m.cnf"]);
  status("Compiling...");
  const r2 = await D.runTool("ganak", ["--verb", "1", "--mode", "1", "--compile", "/out.bdnnf", "/m.cnf"],
    { "/m.cnf": cnf }, ["/out.bdnnf"]);
  if (g !== S.gen) return;
  D.log("ganak --mode 1 --compile out.bdnnf m.cnf", r2);
  const cm = r2.stdout.match(/c s exact arb \S+ (\S+)/);
  if (r2.status !== 0 || !r2.files["/out.bdnnf"] || !cm) {
    showErrors([], "Ganak failed:\n" + r2.stderr.split("\n").slice(-5).join("\n"));
    status("");
    $("gm-results").classList.add("stale");
    return;
  }
  // Variable numbers change with the model, so the model, its CNF and its
  // circuit are only ever replaced together.
  S.model = model;
  S.cnf = cnf;
  S.circuit = r2.files["/out.bdnnf"];
  S.prior = null;
  S.baseLits = null;
  for (const k of [...S.evidence.keys()]) if (!model.nodes.some((n) => n.name === k)) S.evidence.delete(k);
  S.z = Q.parse(cm[1]);
  const compileMs = r2.ms;
  status("Answering...");
  const r3 = await bdnnf(["batch", "/c.bdnnf", "/q.txt", "--cnf", "/m.cnf"], { "/q.txt": batchFile(model) });
  const r4 = await bdnnf(["lits", "/c.bdnnf", "--cnf", "/m.cnf"]);
  if (g !== S.gen) return;
  S.prior = margMap(parseLits(r4.stdout)).marg;
  renderAnswers(answers(model, r3.stdout), compileMs);
  fillTops();
  await updateEvidence(g);
  $("gm-results").classList.remove("stale");
  status("");
}

function renderAnswers(ans: Answer[], compileMs: number) {
  const m = S.model!;
  const box = $("gm-answers");
  box.innerHTML = "";
  const stats = el("div", { class: "stats" },
    stat("events", String(randomNodes(m).filter((n) => !n.anon).length)),
    stat("definitions", String(m.nodes.filter((n) => n.kind === "def").length)),
    stat("compile", `${compileMs.toFixed(0)} ms`));
  box.append(stats);
  if (m.stmts.some((s) => s.kind === "constraint"))
    box.append(el("p", { class: "muted" }, "The constraints hold with probability ", pq(S.z),
      "; every answer is conditioned on them."));
  if (!ans.length) {
    box.append(el("p", { class: "muted" }, "No query statements: add e.g. ", el("code", {}, "query lost"), "."));
    return;
  }
  const rows = ans.map((a) => el("tr", {},
    el("td", {}, el("a", { href: "#", class: "line", "data-line": String(a.stmt.line) },
      pretty(a.stmt) + (a.stmt.kind === "require" ? ` ${a.stmt.op} ${fmtP(Q.parse(a.stmt.bound!))}` : ""))),
    el("td", { class: "num" }, pq(a.p, "big-p")),
    el("td", {}, a.stmt.kind === "require"
      ? el("span", { class: "badge " + (a.pass ? "pass" : "fail") }, a.pass ? "PASS" : "FAIL") : "")));
  const t = el("table", { class: "list answers" },
    el("thead", {}, el("tr", {}, el("th", {}, "question"), el("th", {}, "probability"), el("th", {}, ""))),
    el("tbody", {}, ...rows));
  t.querySelectorAll<HTMLAnchorElement>("a.line").forEach((a) => a.addEventListener("click", (e) => {
    e.preventDefault();
    selectLine(Number(a.dataset.line));
  }));
  box.append(t);
}

function stat(k: string, v: string): HTMLElement {
  return el("div", { class: "stat" }, el("div", { class: "muted" }, k), el("div", { class: "stat-v", title: v }, v));
}

async function updateEvidence(g = S.gen) {
  const e = ++S.ev;
  const ev = evidenceLits();
  const r = await bdnnf(["lits", "/c.bdnnf", "--cnf", "/m.cnf", ...assumeArgs(ev)], {}, ev.length > 0);
  if (g !== S.gen || e !== S.ev) return;
  S.baseLits = parseLits(r.stdout);
  S.baseEv = e;
  const { marg, total } = margMap(S.baseLits);
  S.pe = S.z.isZero() ? null : total.div(S.z);
  S.marg = marg;
  const box = $("gm-evidence");
  box.innerHTML = "";
  if (!S.evidence.size) box.append(el("span", { class: "muted" }, "No evidence: click a node to assume it true, again for false, once more to forget it."));
  else {
    for (const [name, v] of S.evidence) {
      const chip = el("button", { class: "chip " + (v ? "core" : "dead"), title: "remove" }, `${v ? "" : "¬"}${name} ✕`);
      chip.addEventListener("click", () => {
        S.evidence.delete(name);
        updateEvidence();
      });
      box.append(chip);
    }
    const clear = el("button", {}, "clear");
    clear.addEventListener("click", () => {
      S.evidence.clear();
      updateEvidence();
    });
    box.append(clear, el("span", { class: "muted" }, " P(evidence) = "),
      S.pe && !S.pe.isZero() ? pq(S.pe) : el("span", { class: "error" }, "0: impossible evidence"));
  }
  if (!await drawGraph(e)) return;
  await whatMatters(g);
}

async function drawGraph(e: number): Promise<boolean> {
  const m = S.model!;
  const zeroEv = S.evidence.size > 0 && (!S.pe || S.pe.isZero());
  const dot = modelDot(m, zeroEv ? null : S.marg, S.prior, S.evidence);
  const g = await D.renderDot(dot);
  if (e !== S.ev) return false;
  g.querySelectorAll<SVGGElement>("g.node").forEach((node) => {
    const i = Number((node.id || "").replace(/^gm/, ""));
    const n = m.nodes[i];
    if (!n) return;
    node.style.cursor = "pointer";
    node.addEventListener("click", () => {
      const v = S.evidence.get(n.name);
      if (v === undefined) S.evidence.set(n.name, true);
      else if (v) S.evidence.set(n.name, false);
      else S.evidence.delete(n.name);
      updateEvidence();
    });
  });
  const box = $("gm-graph");
  box.innerHTML = "";
  box.append(g);
  return true;
}

function fillTops() {
  const m = S.model!;
  const sel = $<HTMLSelectElement>("gm-top");
  const old = sel.value;
  sel.innerHTML = "";
  const defs = m.nodes.filter((n) => n.kind === "def" && !n.anon);
  const queried = new Set(m.stmts.filter((s) => s.kind !== "constraint").map((s) => s.lit));
  for (const n of defs) sel.append(el("option", { value: n.name }, n.label ? `${n.label} (${n.name})` : n.name));
  const def = defs.find((n) => n.name === old) ?? defs.find((n) => queried.has(n.var)) ?? defs[defs.length - 1];
  if (def) sel.value = def.name;
}

async function whatMatters(g = S.gen) {
  const box = $("gm-analysis");
  const top = S.model && S.baseLits && S.baseEv === S.ev ? nodeByName($<HTMLSelectElement>("gm-top").value) : undefined;
  if (!top) {
    box.innerHTML = "";
    return;
  }
  const e = S.ev, w = ++S.wm;
  const r = await bdnnf(["lits", "/c.bdnnf", "--cnf", "/m.cnf", ...assumeArgs([...evidenceLits(), top.var])], {}, false);
  // the base marginals must be for the same evidence as this run
  if (g !== S.gen || e !== S.ev || w !== S.wm || S.baseEv !== e) return;
  const { pT, rows } = impacts(S.model!, S.baseLits!, parseLits(r.stdout));
  const name = `“${top.label || top.name}”`;
  box.innerHTML = "";
  if (!pT) {
    box.append(el("p", { class: "error" }, "The evidence is impossible."));
    return;
  }
  box.append(el("p", {}, `P(${top.name}) = `, pq(pT, "big-p"), S.evidence.size ? el("span", { class: "muted" }, " given the evidence") : ""));
  const cell = (q: Q | null, cls: string) => q ? el("td", { class: "num " + cls }, pq(q)) : el("td", { class: "muted" }, "impossible");
  const up = (q: Q | null) => q && !pT.isZero() && q.toNumber() >= 10 * pT.toNumber() ? "hi-up" : "";
  const down = (q: Q | null) => q && (pT.isZero() ? false : q.toNumber() <= pT.toNumber() / 10) ? "hi-down" : "";
  box.append(el("div", { class: "scroll" }, el("table", { class: "list matters" },
    el("thead", {}, el("tr", {}, el("th", {}, "event"), el("th", {}, "how likely"),
      el("th", {}, `P(${top.name}) if it happens`), el("th", {}, `P(${top.name}) if it never happens`))),
    el("tbody", {}, ...rows.map((x) => el("tr", {},
      el("td", { title: x.node.name }, x.node.label || x.node.name),
      el("td", { class: "num" }, pq(x.p)), cell(x.ifYes, up(x.ifYes)), cell(x.ifNo, down(x.ifNo))))))));
  box.append(el("p", { class: "muted" }, `The events that move ${name} most come first. Red: it becomes at least 10 times likelier; `
    + "green: at least 10 times less likely. Everything follows the evidence set in the graph."));
}

function setText(t: string) {
  $<HTMLTextAreaElement>("gm-text").value = t;
  S.evidence.clear();
  $("gm-analysis").innerHTML = "";
  compute();
}

export function initModelTab(deps: Deps) {
  D = deps;
  const sel = $<HTMLSelectElement>("gm-examples");
  MODEL_EXAMPLES.forEach((x, i) => sel.append(el("option", { value: String(i) }, x.name)));
  sel.addEventListener("change", () => setText(MODEL_EXAMPLES[Number(sel.value)].text));
  const ta = $<HTMLTextAreaElement>("gm-text");
  ta.addEventListener("input", () => {
    $("gm-results").classList.add("stale");
    clearTimeout(timer);
    timer = window.setTimeout(compute, 700);
  });
  ta.addEventListener("keydown", (e) => {
    if (e.key === "Tab") {
      e.preventDefault();
      ta.setRangeText("  ", ta.selectionStart, ta.selectionEnd, "end");
    }
  });
  $<HTMLInputElement>("gm-file").addEventListener("change", async (e) => {
    const f = (e.target as HTMLInputElement).files?.[0];
    if (f) setText(await f.text());
  });
  $("gm-compute").addEventListener("click", compute);
  $("gm-cnf").addEventListener("click", () => {
    if (S.cnf) D.openInCnf(S.cnf);
  });
  $("gm-download").addEventListener("click", () => {
    const a = el("a", { href: URL.createObjectURL(new Blob([$<HTMLTextAreaElement>("gm-text").value])), download: "model.gm" });
    a.click();
  });
  $("gm-top").addEventListener("change", () => whatMatters());
  setText(MODEL_EXAMPLES[0].text);
}
