#!/usr/bin/env python3
# Fuzz the bounded table counters (--frontierstates, --groupedsweep) against
# normal counting. Random CNFs shaped to reach each method: small random
# k-CNFs and binary-heavy residuals (word table, bucket elimination, clause
# frontier), banded equal-scope clause groups (grouped root and residual
# sweep), and selector blocks (grouped root selector count). Random budgets,
# with and without Arjun, sometimes projected (the counters must decline).
# Usage: scripts/fuzz_sweep.py [--num N] [--ganak build/ganak] [-j 16]
import argparse, os, random, re, subprocess, sys
from multiprocessing.pool import ThreadPool as Pool  # workers only wait on subprocesses

def kcnf(rng, base, n):
    cls = []
    ratio = rng.uniform(1.0, 4.5)
    for _ in range(int(ratio * n)):
        k = rng.choice([2, 2, 3, 3, 3, 4]) if rng.random() < 0.5 else rng.choice([2, 3])
        vs = rng.sample(range(n), min(k, n))
        cls.append([(base + v + 1) * rng.choice([-1, 1]) for v in vs])
    return cls

def banded(rng, base, n):
    """Clauses grouped on shared scopes along a band: small frontier, many groups."""
    cls = []
    width = rng.randint(2, 5)
    for start in range(0, max(1, n - width + 1), rng.randint(1, width)):
        scope = list(range(start, min(n, start + width)))
        if rng.random() < 0.3 and start + width + 2 < n:  # occasional long edge
            scope[-1] = rng.randrange(start + width, n)
        for _ in range(rng.randint(1, 2 ** len(scope) - 1)):
            cls.append([(base + v + 1) * rng.choice([-1, 1]) for v in scope])
    return cls

def selectors(rng, base):
    """Selector blocks as the grouped root recognizes them: blocks of equally
    many positive pairs; each selector picks one side of every pair of its
    block; each row clause holds one selector per block; negative clauses
    forbid one side-choice of a whole block. Sometimes slightly off-shape."""
    cls, v = [], base
    nblocks, width = rng.randint(2, 4), rng.randint(1, 3)
    blocks = []
    for _ in range(nblocks):
        pairs = []
        for _ in range(width):
            cls.append([v + 1, v + 2]); pairs.append((v + 1, v + 2)); v += 2
        blocks.append(pairs)
    for _ in range(rng.randint(1, 4)):
        row = []
        for pairs in blocks:
            v += 1
            for pr in pairs: cls.append([rng.choice(pr), -v])
            row.append(v)
        cls.append(row)
    if width >= 2:
        for _ in range(rng.randint(0, 2 * nblocks)):
            cls.append([-rng.choice(pr) for pr in rng.choice(blocks)])
    if rng.random() < 0.2:  # off-shape: must decline or still agree
        a, b = rng.sample(range(base + 1, v + 1), 2)
        cls.append([rng.choice([-1, 1]) * a, rng.choice([-1, 1]) * b])
    return cls, v - base

def build(rng):
    cls, nv = [], 0
    if rng.random() < 0.2:
        cls, nv = selectors(rng, 0)
        return nv, cls
    for _ in range(rng.choice([1, 1, 2, 3])):  # disconnected parts give several components
        kind = rng.randrange(3)
        if kind == 0:
            n = rng.randint(6, 60); cls += kcnf(rng, nv, n)
        elif kind == 1:
            n = rng.randint(8, 70); cls += banded(rng, nv, n)
        else:  # binary-heavy with a few long clauses
            n = rng.randint(20, 90)
            c = kcnf(rng, nv, n)
            cls += [x[:2] if len(x) > 2 and rng.random() < 0.7 else x for x in c]
        nv += n
    if rng.random() < 0.15 and nv > 1:  # some free variables
        nv += rng.randint(1, 3)
    return max(nv, 1), cls

def run(ganak, opts, fname, timeout):
    try:
        out = subprocess.run([ganak] + opts + [fname], capture_output=True, text=True, timeout=timeout).stdout
    except subprocess.TimeoutExpired:
        return None, ""
    m = re.search(r"^c s exact arb int (\d+)", out, re.M)
    return (int(m.group(1)) if m else "NOCOUNT"), out

def one(args):
    seed, ganak, tmpdir, timeout = args
    rng = random.Random(seed)
    nv, cls = build(rng)
    rng.shuffle(cls)
    fname = os.path.join(tmpdir, f"sweep-{seed}.cnf")
    with open(fname, "w") as f:
        f.write(f"p cnf {nv} {len(cls)}\n")
        if rng.random() < 0.15:
            show = sorted(rng.sample(range(1, nv + 1), rng.randint(1, nv)))
            f.write("c p show " + " ".join(map(str, show)) + " 0\n")
        for c in cls: f.write(" ".join(map(str, c)) + " 0\n")
    common = ["--arjun", str(rng.choice([0, 0, 1])), "--td", str(rng.choice([0, 1]))]
    opts = common + [
        "--frontierstates", str(rng.choice([0, 1, 2, 64, 8192, 8192, 1 << 20])),
        "--frontierminvars", str(rng.choice([4, 4, 12, 30])),
        "--frontierbucket", str(rng.choice([0, 1, 1])),
        "--groupedsweep", str(rng.choice([0, 1]))]
    got, out = run(ganak, opts, fname, timeout)
    ref, _ = run(ganak, common + ["--frontierstates", "0", "--groupedsweep", "0"], fname, timeout)
    msg = None
    if got == "NOCOUNT" or ref == "NOCOUNT":
        msg = f"no count / crash: {out[-500:]}"
    elif got is not None and ref is not None and got != ref:
        msg = f"WRONG: sweep {got} vs normal {ref}"
    used = set()
    if "grouped root selector" in out: used.add("root-selector")
    if "grouped root scopes" in out: used.add("root-grouped")
    m = re.search(r"residual attempts \d+ finished \d+ \(words (\d+) buckets (\d+) sparse (\d+)", out)
    if m:
        for name, k in zip(("word", "bucket", "sparse"), m.groups()):
            if int(k): used.add(name)
    if re.search(r"grouped residual attempts \d+ finished [1-9]", out): used.add("residual-grouped")
    if msg is None:
        os.unlink(fname)
        return seed, None, used
    return seed, f"{msg}\n  ganak {' '.join(opts)} {fname}", used

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--num", type=int, default=300)
    p.add_argument("--ganak", default="build/ganak")
    p.add_argument("--seed", type=int, default=random.randrange(1 << 30))
    p.add_argument("-j", type=int, default=16)
    p.add_argument("-t", type=float, default=30)
    a = p.parse_args()
    tmpdir = "fuzz-sweep-out"
    os.makedirs(tmpdir, exist_ok=True)
    fails = used = 0
    paths = {}
    with Pool(a.j) as pool:
        for seed, err, u in pool.imap_unordered(one, [(a.seed+i, a.ganak, tmpdir, a.t) for i in range(a.num)]):
            used += bool(u)
            for k in u: paths[k] = paths.get(k, 0) + 1
            if err:
                fails += 1
                print(f"seed {seed}: {err}", flush=True)
    print(f"{a.num} runs, {used} used a table counter {dict(sorted(paths.items()))}, {fails} failures (base seed {a.seed})")
    sys.exit(1 if fails else 0)

if __name__ == "__main__":
    main()
