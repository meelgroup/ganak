#!/usr/bin/env python3
# Fuzz ganak-model against brute-force enumeration. Usage: scripts/fuzz_gmodel.py [--num N] [-j J] [--seed S]

import argparse
import itertools
import json
import os
import random
import re
import subprocess
import sys
import tempfile
import traceback
from fractions import Fraction as F
from multiprocessing import Pool

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BUILD = os.path.join(ROOT, "build")
GM = os.path.join(BUILD, "ganak-model")
GANAK = os.path.join(BUILD, "ganak")
BDNNF = os.path.join(BUILD, "bdnnf")
OUT = os.path.join(ROOT, "fuzz-gmodel-out")
MAX_RANDOM = 12


class Fail(Exception):
    pass


def run(cmd, inp=None, timeout=60):
    r = subprocess.run(cmd, input=inp, capture_output=True, text=True, timeout=timeout)
    if r.returncode < 0 or r.returncode in (134, 136, 139):
        raise Fail(f"signal/crash {r.returncode}: {' '.join(cmd)}\n{r.stderr[-2000:]}")
    return r


def prob_text(rnd, p):
    if p in (0, 1):
        return rnd.choice([str(p), f"{p}.0", f"{p * 100}%"])
    if (p * 100).denominator == 1 and rnd.random() < 0.3:
        return f"{p * 100}%"
    if 10 ** 6 % p.denominator == 0 and rnd.random() < 0.6:
        d = p.numerator * (10 ** 6 // p.denominator)
        s = f"{d / 10**6:.6f}".rstrip("0")
        if rnd.random() < 0.4:
            m, e = d, -6
            while m % 10 == 0:
                m //= 10
                e += 1
            return f"{m}e{e}"
        return s[1:] if s.startswith("0.") and rnd.random() < 0.3 else s
    return f"{p.numerator}/{p.denominator}"


def rand_prob(rnd):
    r = rnd.random()
    if r < 0.05:
        return F(0)
    if r < 0.1:
        return F(1)
    return rnd.choice([F(1, 2), F(1, 4), F(3, 4), F(1, 10), F(1, 5), F(2, 5), F(1, 8), F(1, 1000), F(1, 3), F(2, 7)])


PREC = {"ite": 0, "iff": 1, "imp": 2, "or": 3, "xor": 4, "and": 5, "not": 6}


class Model:
    def __init__(self, rnd):
        self.rnd = rnd
        self.random = []  # (key, p, name or None if anonymous) ; key used in worlds
        self.choices = []  # (name, [(value, p)])
        self.items = []  # defs in dependency order: (name, kind, data)
        self.events = []
        self.stmts = []
        self.flipn = 0

    def new_flip(self, p):
        self.flipn += 1
        key = f"#f{self.flipn}"
        self.random.append((key, p, None))
        return ("flip", key, p)

    def atoms(self, upto):
        a = [("ref", e) for e in self.events]
        for c, vals in self.choices:
            a += [("ref", f"{c}.{v}") for v, _ in vals]
        a += [("ref", it[0]) for it in self.items[:upto]]
        return a

    def expr(self, upto, depth):
        rnd = self.rnd
        atoms = self.atoms(upto)
        if depth <= 0 or rnd.random() < 0.25:
            r = rnd.random()
            if r < 0.05 or not atoms:
                return ("const", rnd.random() < 0.5)
            if r < 0.12 and len(self.random) < MAX_RANDOM:
                return self.new_flip(rand_prob(rnd))
            return rnd.choice(atoms)
        op = rnd.choice(["not", "and", "or", "xor", "iff", "imp", "ite", "card", "and", "or", "not"])
        if op == "not":
            return ("not", self.expr(upto, depth - 1))
        if op == "ite":
            return ("ite", self.expr(upto, depth - 1), self.expr(upto, depth - 1), self.expr(upto, depth - 1))
        if op == "card":
            n = rnd.randint(1, 4)
            kids = [self.expr(upto, depth - 2) for _ in range(n)]
            return ("card", rnd.choice(["atleast", "atmost", "exactly"]), rnd.randint(0, n + 1), kids)
        return (op, self.expr(upto, depth - 1), self.expr(upto, depth - 1))

    def show(self, e, ctx=0):
        rnd = self.rnd
        t = e[0]
        if t == "const":
            return "true" if e[1] else "false"
        if t == "ref":
            return e[1]
        if t == "flip":
            return f"flip({prob_text(rnd, e[2])})"
        if t == "card":
            of = " of" if rnd.random() < 0.5 else ""
            s = f"{e[1]} {e[2]}{of} (" + ", ".join(self.show(k, 0) for k in e[3]) + ")"
            return s
        lvl = PREC[t]
        if t == "not":
            s = rnd.choice(["!", "~", "not "]) + self.show(e[1], 6)
        elif t == "ite":
            s = f"if {self.show(e[1], 0)} then {self.show(e[2], 0)} else {self.show(e[3], 0)}"
        else:
            sym = {"and": ["&", "and"], "or": ["|", "or"], "xor": ["^", "xor"], "iff": ["<->"], "imp": ["->"]}[t]
            sep = rnd.choice(sym)
            nl = "\n    " if rnd.random() < 0.1 else " "
            if t == "imp":
                s = f"{self.show(e[1], lvl + 1)} {sep}{nl}{self.show(e[2], lvl)}"
            else:
                s = f"{self.show(e[1], lvl)} {sep}{nl}{self.show(e[2], lvl + 1)}"
        if lvl < ctx or (lvl == 0 and ctx > 0) or rnd.random() < 0.08:
            s = "(" + s + ")"
        return s


def gen_model(rnd):
    m = Model(rnd)
    nev = rnd.randint(1, 5)
    for i in range(nev):
        name = f"e{i}"
        p = rand_prob(rnd)
        m.events.append(name)
        m.random.append((name, p, name))
    if rnd.random() < 0.3:
        vals = []
        n = rnd.randint(2, 3)
        ps = [F(rnd.randint(0, 4)) for _ in range(n)]
        if sum(ps) == 0:
            ps[0] = F(1)
        ps = [x / sum(ps) for x in ps]
        for k in range(n):
            vals.append((f"v{k}", ps[k]))
        m.choices.append(("c0", vals))
    nitems = rnd.randint(1, 6)
    for j in range(nitems):
        if rnd.random() < 0.2:
            leak = rnd.random() < 0.3
            name = f"l{j}" if leak else f"r{j}"
            rules = []
            for _ in range(rnd.randint(1, 3)):
                p = rand_prob(rnd) if rnd.random() < 0.6 else F(1)
                body = m.expr(j, 2)
                key = None
                if p != 1:
                    m.flipn += 1
                    key = f"#r{m.flipn}"
                    m.random.append((key, p, None))
                rules.append((body, p, key))
            lp = None
            if leak:
                lp = rand_prob(rnd)
                m.random.append((name + ".self", lp, name + ".self"))
            m.items.append((name, "rules", (rules, lp)))
        else:
            m.items.append((f"d{j}", "def", m.expr(j, rnd.randint(1, 4))))
    allx = len(m.items)
    for _ in range(rnd.randint(0, 2) if rnd.random() < 0.4 else 0):
        m.stmts.append(("constraint", m.expr(allx, 2)))
    for _ in range(rnd.randint(1, 4)):
        given = [m.expr(allx, 1) for _ in range(rnd.randint(0, 2) if rnd.random() < 0.5 else 0)]
        m.stmts.append(("query", m.expr(allx, 3), given))
    if rnd.random() < 0.4:
        given = [m.expr(allx, 1)] if rnd.random() < 0.3 else []
        m.stmts.append(("require", m.expr(allx, 2), given, rnd.choice(["<", "<=", ">", ">="]), rand_prob(rnd)))
    return m


def model_text(m):
    rnd = m.rnd
    lines = []
    leak_p = {}
    for name, kind, data in m.items:
        if kind == "rules" and data[1] is not None:
            leak_p[name] = data[1]
    decls = []
    for e in m.events:
        p = next(x[1] for x in m.random if x[0] == e)
        lab = f' "Event {e}"' if rnd.random() < 0.3 else ""
        decls.append(f"event {e}{lab} {prob_text(rnd, p)}")
    for name, p in leak_p.items():
        decls.append(f"event {name} {prob_text(rnd, p)}")
    for c, vals in m.choices:
        decls.append(f"choice {c} {{ " + ", ".join(f"{v} {prob_text(rnd, p)}" for v, p in vals) + " }")
    for name, kind, data in m.items:
        if kind == "def":
            lab = f' "Def {name}"' if rnd.random() < 0.3 else ""
            decls.append(f"{name}{lab} = {m.show(data, 0)}")
        else:
            for body, p, key in data[0]:
                ps = "" if p == 1 and rnd.random() < 0.5 else " " + prob_text(rnd, p)
                decls.append(f"rule {m.show(body, 3)} -> {name}{ps}")
    rnd.shuffle(decls)
    lines += decls
    for s in m.stmts:
        if s[0] == "constraint":
            lines.append(f"constraint {m.show(s[1], 0)}")
        else:
            g = " given " + ", ".join(m.show(x, 0) for x in s[2]) if s[2] else ""
            if s[0] == "query":
                lines.append(f"query {m.show(s[1], 0)}{g}")
            else:
                lines.append(f"require {m.show(s[1], 0)}{g} {s[3]} {prob_text(rnd, s[4])}")
    out = []
    for ln in lines:
        if rnd.random() < 0.1:
            out.append("# a comment")
        out.append(ln + (" # trailing" if rnd.random() < 0.1 else ""))
        if rnd.random() < 0.1:
            out.append("")
    return "\n".join(out) + "\n"


def ev(e, val, flips):
    t = e[0]
    if t == "const":
        return e[1]
    if t == "ref":
        return val[e[1]]
    if t == "flip":
        return flips[e[1]]
    if t == "not":
        return not ev(e[1], val, flips)
    if t == "and":
        return ev(e[1], val, flips) and ev(e[2], val, flips)
    if t == "or":
        return ev(e[1], val, flips) or ev(e[2], val, flips)
    if t == "xor":
        return ev(e[1], val, flips) != ev(e[2], val, flips)
    if t == "iff":
        return ev(e[1], val, flips) == ev(e[2], val, flips)
    if t == "imp":
        return (not ev(e[1], val, flips)) or ev(e[2], val, flips)
    if t == "ite":
        return ev(e[2], val, flips) if ev(e[1], val, flips) else ev(e[3], val, flips)
    if t == "card":
        c = sum(1 for k in e[3] if ev(k, val, flips))
        return c >= e[2] if e[1] == "atleast" else c <= e[2] if e[1] == "atmost" else c == e[2]
    raise ValueError(t)


def worlds(m):
    bools = [(k, p) for k, p, _ in m.random]
    for bits in itertools.product([False, True], repeat=len(bools)):
        w = F(1)
        assign = {}
        for (k, p), b in zip(bools, bits):
            w *= p if b else 1 - p
            assign[k] = b
        if w == 0:
            continue
        for cvals in itertools.product(*[vals for _, vals in m.choices]):
            ww = w
            for (c, vals), (v, p) in zip(m.choices, cvals):
                ww *= p
                for vv, _ in vals:
                    assign[f"{c}.{vv}"] = vv == v
            if ww == 0:
                continue
            yield ww, dict(assign)


def values(m, assign):
    val = dict(assign)
    for name, kind, data in m.items:
        if kind == "def":
            val[name] = ev(data, val, assign)
        else:
            rules, leak = data
            x = bool(leak is not None and assign[name + ".self"])
            for body, p, key in rules:
                if ev(body, val, assign) and (key is None or assign[key]):
                    x = True
            val[name] = x
    return val


def expected(m):
    z = F(0)
    stmts = [s for s in m.stmts if s[0] != "constraint"]
    num = [F(0)] * len(stmts)
    den = [F(0)] * len(stmts)
    marg = {}
    for w, a in worlds(m):
        val = values(m, a)
        if not all(ev(s[1], val, a) for s in m.stmts if s[0] == "constraint"):
            continue
        z += w
        for k, v in val.items():
            if v and not k.startswith("#"):
                marg[k] = marg.get(k, F(0)) + w
        for i, s in enumerate(stmts):
            if all(ev(g, val, a) for g in s[2]):
                den[i] += w
                if ev(s[1], val, a):
                    num[i] += w
    return z, num, den, marg


def parse_solve(out):
    res = []
    for ln in out.splitlines():
        if ln.startswith("P(constraints)"):
            continue
        if " is undefined" in ln:
            res.append(None)
            continue
        mm = re.search(r"= (\S+) ~", ln) if ln.startswith("P(") else re.search(r": (PASS|FAIL) \((\S+) ~", ln)
        if mm:
            res.append(F(mm.group(mm.lastindex)))
    return res


def check_model(rnd, tmp, idx):
    m = gen_model(rnd)
    text = model_text(m)
    path = os.path.join(tmp, f"m{idx}.gm")
    with open(path, "w") as f:
        f.write(text)
    z, num, den, marg = expected(m)
    stmts = [s for s in m.stmts if s[0] != "constraint"]
    r = run([GM, "--solve", path])
    if r.returncode == 1 and z != 0:
        raise Fail(f"error exit:\n{r.stderr}{r.stdout}", ) if r.stderr else Fail("exit 1:\n" + r.stdout)
    if z == 0:
        if r.returncode != 1 or "never hold" not in r.stdout:
            raise Fail(f"Z=0 expected 'never hold', got exit {r.returncode}:\n{r.stdout}{r.stderr}")
        return path
    got = parse_solve(r.stdout)
    if len(got) != len(stmts):
        raise Fail(f"expected {len(stmts)} answers, got:\n{r.stdout}{r.stderr}")
    fails = False
    for i, s in enumerate(stmts):
        exp = None if den[i] == 0 else num[i] / den[i]
        if got[i] != exp:
            raise Fail(f"statement {i + 1}: expected {exp}, got {got[i]}\n{r.stdout}")
        if s[0] == "require" and exp is not None:
            op, b = s[3], s[4]
            ok = exp < b if op == "<" else exp <= b if op == "<=" else exp > b if op == ">" else exp >= b
            fails |= not ok
            answer_lines = [ln for ln in r.stdout.splitlines() if ln.startswith(("P(", "require "))
                            and not ln.startswith("P(constraints)")]
            if ("PASS" if ok else "FAIL") not in answer_lines[i]:
                raise Fail(f"requirement {i + 1} verdict wrong\n{r.stdout}")
        if s[0] == "require" and exp is None:
            fails = True
    want = 3 if fails else 0
    if r.returncode != want:
        raise Fail(f"exit {r.returncode}, expected {want}\n{r.stdout}")

    if rnd.random() < 0.5:
        check_compiled(m, path, tmp, idx, z, num, den, marg, stmts)
    return path


def check_compiled(m, path, tmp, idx, z, num, den, marg, stmts):
    cnf = os.path.join(tmp, f"m{idx}.cnf")
    js = os.path.join(tmp, f"m{idx}.json")
    run([GM, path, "-o", cnf, "--json", js])
    with open(js) as f:
        model = json.load(f)
    circ = os.path.join(tmp, f"m{idx}.bdnnf")
    r = run([GANAK, "--mode", "1", "--compile", circ, cnf])
    mm = re.search(r"c s exact arb \S+ (\S+)", r.stdout)
    if not mm or F(mm.group(1)) != z:
        raise Fail(f"ganak --compile count {mm and mm.group(1)} != {z}\n{r.stdout[-1500:]}")
    lits = run([BDNNF, "lits", circ, "--cnf", cnf])
    cnt = {}
    for ln in lits.stdout.splitlines():
        if ln.startswith("c") or not ln.strip():
            continue
        a, b = ln.split()
        cnt[int(a)] = F(b)
    for n in model["nodes"]:
        if n["anon"] and not n["name"].endswith(".self"):
            continue
        v = n["var"]
        if cnt.get(v, F(0)) != marg.get(n["name"], F(0)):
            raise Fail(f"bdnnf lits: {n['name']} (var {v}) {cnt.get(v)} != {marg.get(n['name'], 0)}")
    js_stmts = [s for s in model["stmts"] if s["kind"] != "constraint"]
    for i, s in enumerate(js_stmts):
        g = s["given"]
        def count(ls):
            args = [BDNNF, "count", circ, "--cnf", cnf]
            if ls:
                args += ["--assume", ",".join(str(x) for x in ls)]
            out = run(args).stdout.split()
            return F(out[-1]) if out else F(0)
        d = count(g)
        if d != den[i]:
            raise Fail(f"bdnnf count given of statement {i + 1}: {d} != {den[i]}")
        if d != 0 and count(g + [s["lit"]]) != num[i]:
            raise Fail(f"bdnnf count of statement {i + 1} != {num[i]}")


def check_mutation(rnd, tmp, idx):
    m = gen_model(rnd)
    text = list(model_text(m))
    for _ in range(rnd.randint(1, 5)):
        op = rnd.random()
        pos = rnd.randrange(len(text) + 1)
        if op < 0.4 and text:
            del text[min(pos, len(text) - 1)]
        elif op < 0.8:
            text.insert(pos, rnd.choice('()!&|^-><=,;{}"#\n .01e%/abq' + "\t"))
        else:
            text.insert(pos, rnd.choice([" event ", " rule ", " given ", " if ", " then ", " atleast ", "flip(", " query "]))
    src = "".join(text)
    path = os.path.join(tmp, f"x{idx}.gm")
    with open(path, "w") as f:
        f.write(src)
    r = run([GM, "--check", "--json", "-", path])
    if r.returncode not in (0, 1):
        raise Fail(f"mutation: exit {r.returncode}\n{r.stderr}")
    js = json.loads(r.stdout)
    if (r.returncode == 1) != bool(js["errors"]):
        raise Fail(f"mutation: exit {r.returncode} but errors {js['errors']}")
    if r.returncode == 0 and rnd.random() < 0.5:
        run([GM, "--solve", path])
    return path


def one(args):
    seed, kinds = args
    rnd = random.Random(seed)
    kind = rnd.choice(kinds)
    with tempfile.TemporaryDirectory() as tmp:
        try:
            {"model": check_model, "mutate": check_mutation}[kind](
                rnd, tmp, seed)
            return None
        except Exception as e:
            os.makedirs(OUT, exist_ok=True)
            for fn in os.listdir(tmp):
                if fn.endswith((".gm", ".pl", ".xml")):
                    with open(os.path.join(tmp, fn)) as src, open(os.path.join(OUT, f"{seed}-{fn}"), "w") as dst:
                        dst.write(src.read())
            msg = str(e) if isinstance(e, Fail) else traceback.format_exc()
            return f"seed {seed} ({kind}): {msg}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--num", type=int, default=500)
    ap.add_argument("-j", type=int, default=16)
    ap.add_argument("--seed", type=int, default=None)
    ap.add_argument("--kinds", default="model,model,model,mutate")
    a = ap.parse_args()
    base = a.seed if a.seed is not None else random.randrange(10 ** 9)
    kinds = a.kinds.split(",")
    jobs = [(base + i, kinds) for i in range(a.num)]
    bad = 0
    with Pool(a.j) as pool:
        for i, res in enumerate(pool.imap_unordered(one, jobs)):
            if res:
                bad += 1
                print(res, flush=True)
    print(f"{a.num} tests from seed {base}, {bad} failures" + (f" (files in {OUT})" if bad else ""))
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
