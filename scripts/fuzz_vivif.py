#!/usr/bin/env python3
# Stress-fuzz in-search vivification: build-vivif/ganak (built with -DCHECK_VIVIF) vs build/ganak --vivif 0.
# Usage: scripts/fuzz_vivif.py [--num N] [-j J] [--tout T] [--seed S]
import argparse
import os
import random
import re
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BIERE = os.path.join(ROOT, "..", "count_fuzzer", "generators", "biere-fuzz")
CHK = os.path.join(ROOT, "build-vivif", "ganak")
REF = os.path.join(ROOT, "build", "ganak")
OUT = os.path.join(ROOT, "fuzz-vivif-out")


def gen_random_cnf(rnd, path, proj):
    # Near-threshold mixed-width CNF with a few blocks, so the search hits conflicts
    # deep in the tree and components split.
    n = rnd.randint(25, 140)
    blocks = rnd.choice([1, 1, 2, 3, 5])
    cls = []
    for b in range(blocks):
        lo = 1 + b * n // blocks
        hi = (b + 1) * n // blocks
        width = hi - lo + 1
        if width < 6:
            continue
        ratio = rnd.uniform(2.5, 4.6)
        for _ in range(int(width * ratio)):
            k = rnd.choice([2, 3, 3, 3, 4, 4, 5, 7])
            vs = rnd.sample(range(lo, hi + 1), min(k, width))
            cls.append([v if rnd.random() < 0.5 else -v for v in vs])
    for _ in range(rnd.randint(0, n // 4)):
        vs = rnd.sample(range(1, n + 1), 3)
        cls.append([v if rnd.random() < 0.5 else -v for v in vs])
    with open(path, "w") as f:
        f.write(f"p cnf {n} {len(cls)}\n")
        if proj:
            show = sorted(rnd.sample(range(1, n + 1), rnd.randint(1, n)))
            f.write("c p show " + " ".join(map(str, show)) + " 0\n")
        for c in cls:
            f.write(" ".join(map(str, c)) + " 0\n")


def stress_opts(rnd):
    o = {
        "vivifevery": rnd.choice(["1", "1", "2", "5", "30"]),
        "vivifsubs": rnd.choice(["0", "1", "1"]),
        "vivifmult": rnd.choice(["0.001", "0.05", "1", "10"]),
        "totusedcutoffvivif": rnd.choice(["0", "5", "50"]),
        "vivifoutern": rnd.choice(["1", "1", "2"]),
        "restart": rnd.choice(["0", "1", "1"]),
        "rstfirst": rnd.choice(["3", "10", "100"]),
        "arjun": rnd.choice(["0", "0", "1"]),
        "td": rnd.choice(["0", "1"]),
        "cache": rnd.choice(["0", "1", "1"]),
        "lbd": rnd.choice(["1", "2", "3"]),
        "rdbclstarget": rnd.choice(["50", "500", "10000"]),
        "rdbeveryn": rnd.choice(["100", "1000", "10000"]),
        "consolidateeveryn": rnd.choice(["100", "1000", "30000"]),
    }
    return [x for k, v in o.items() for x in ("--" + k, v)]


def get_count(out):
    m = re.findall(r"^c s exact \S+ \S+ (.*)$", out, re.M)
    if m:
        return m[-1].strip()
    if re.search(r"^s UNSATISFIABLE", out, re.M):
        return "0"
    return None


def run(cmd, tout):
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=tout)
        return p.returncode, p.stdout + p.stderr
    except subprocess.TimeoutExpired:
        return None, ""


def one(i, seed, tout):
    rnd = random.Random(seed)
    proj = rnd.random() < 0.4
    weighted = rnd.random() < 0.2
    cnf = os.path.join(OUT, f"t{i}_{seed}.cnf")
    if rnd.random() < 0.5 and not weighted:
        gen_random_cnf(rnd, cnf, proj)
    else:
        typ = (1 if proj else 0) + (2 if weighted else 0)
        with open(cnf, "w") as f:
            subprocess.run([BIERE, str(rnd.randint(0, 10**9)), str(typ)], stdout=f, check=True)
    mode = ["--mode", "1"] if weighted else []
    ref_cmd = [REF, "--verb", "0", "--vivif", "0"] + mode + [cnf]
    chk_cmd = [CHK, "--verb", "0"] + mode + stress_opts(rnd) + [cnf]
    rc, out = run(chk_cmd, tout)
    if rc is None:
        os.unlink(cnf)
        return "timeout", None
    bad = None
    if rc != 0 or "ERROR" in out:
        bad = "CRASH/CHECK rc=%s" % rc
    else:
        rrc, rout = run(ref_cmd, tout)
        if rrc is None:
            os.unlink(cnf)
            return "timeout", None
        a, b = get_count(out), get_count(rout)
        if a is None or a != b:
            bad = f"COUNT got {a} expected {b}"
    if bad is None:
        os.unlink(cnf)
        return "ok", None
    tail = "\n".join(l for l in out.splitlines() if not l.startswith("c o"))[-3000:]
    with open(cnf + ".fail", "w") as f:
        f.write(" ".join(chk_cmd) + "\n" + bad + "\n" + tail + "\n")
    return "fail", f"{bad}\n  {' '.join(chk_cmd)}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--num", type=int, default=500)
    ap.add_argument("-j", type=int, default=16)
    ap.add_argument("--tout", type=int, default=30)
    ap.add_argument("--seed", type=int, default=random.randint(0, 10**9))
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    print(f"seed {a.seed}")
    rnd = random.Random(a.seed)
    seeds = [rnd.randint(0, 10**9) for _ in range(a.num)]
    stats = {"ok": 0, "timeout": 0, "fail": 0}
    with ThreadPoolExecutor(a.j) as ex:
        for res, msg in ex.map(lambda x: one(x[0], x[1], a.tout), enumerate(seeds)):
            stats[res] += 1
            if msg:
                print("FAIL:", msg, flush=True)
    print(stats)
    sys.exit(1 if stats["fail"] else 0)


if __name__ == "__main__":
    main()
