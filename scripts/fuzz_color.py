#!/usr/bin/env python3
# Fuzz the direct graph-coloring counter (--colorcnt) against normal counting.
# Usage: scripts/fuzz_color.py [--num N] [--ganak build/ganak] [-j 16]
import argparse, os, random, re, subprocess, sys
from multiprocessing import Pool

def rand_graph(rng, m):
    kind = rng.randrange(4)
    edges = set()
    if kind == 0:
        p = rng.uniform(0.05, 0.6)
        edges = {(u, v) for v in range(m) for u in range(v) if rng.random() < p}
    elif kind == 1:  # sparse, possibly large and disconnected
        for v in range(1, m):
            if rng.random() < 0.85: edges.add((rng.randrange(v), v))
        for _ in range(rng.randrange(m // 4 + 1)):
            u, v = rng.sample(range(m), 2)
            edges.add((min(u, v), max(u, v)))
    elif kind == 2:  # cycle/grid-ish
        for v in range(m):
            edges.add((min(v, (v+1) % m), max(v, (v+1) % m)))
            if v + 5 < m and rng.random() < 0.5: edges.add((v, v+5))
        edges = {e for e in edges if e[0] != e[1]}
    else:  # dense cliquey
        p = rng.uniform(0.6, 1.0)
        edges = {(u, v) for v in range(m) for u in range(v) if rng.random() < p}
    return sorted(edges)

def enc_bits(rng, m, bits, edges):
    k = 1 << bits
    pal = [c for c in range(k) if rng.random() < 0.8] if rng.random() < 0.4 else list(range(k))
    cls = []
    lit = lambda v, j, c: -(v + m*j + 1) if (c >> j) & 1 else (v + m*j + 1)
    for v in range(m):
        for c in range(k):
            if c not in pal: cls.append([lit(v, j, c) for j in range(bits)])
    for u, v in edges:
        for c in range(k):
            if c in pal or rng.random() < 0.3:
                cls.append([lit(w, j, c) for w in (u, v) for j in range(bits)])
    return m*bits, cls

def enc_onehot(m, k, edges):
    var = lambda v, c: v*k + c + 1
    cls = []
    for v in range(m):
        cls.append([var(v, c) for c in range(k)])
        for a in range(k):
            for b in range(a):
                cls.append([-var(v, a), -var(v, b)])
    for u, v in edges:
        for c in range(k):
            cls.append([-var(u, c), -var(v, c)])
    return m*k, cls

def mutate(rng, n, cls):
    r = rng.randrange(5)
    if not cls: r = 1
    if r == 0: del cls[rng.randrange(len(cls))]
    elif r == 1: cls.append([rng.choice([-1, 1]) * rng.randint(1, n)])
    elif r == 2:
        c = cls[rng.randrange(len(cls))]
        c[rng.randrange(len(c))] *= -1
    elif r == 3:
        perm = list(range(1, n+1))
        rng.shuffle(perm)
        cls[:] = [[(1 if l > 0 else -1) * perm[abs(l)-1] for l in c] for c in cls]
    else:
        cls.append(list(cls[rng.randrange(len(cls))]))

def run(ganak, opts, fname, timeout):
    try:
        out = subprocess.run([ganak] + opts + [fname], capture_output=True, text=True, timeout=timeout).stdout
    except subprocess.TimeoutExpired:
        return None, None, ""
    m = re.search(r"^c s exact arb int (\d+)", out, re.M)
    return (int(m.group(1)) if m else "NOCOUNT"), re.search(r"\[color\].*", out), out

def one(args):
    seed, ganak, tmpdir, timeout = args
    rng = random.Random(seed)
    big = rng.random() < 0.15
    m = rng.randint(65, 200) if big else rng.randint(2, 40)
    edges = rand_graph(rng, m)
    if rng.random() < 0.5:
        n, cls = enc_bits(rng, m, rng.choice([2, 2, 3]), edges)
    else:
        n, cls = enc_onehot(m, rng.randint(3, 6), edges)
    mutated = bool(cls) and rng.random() < 0.3
    if mutated:
        for _ in range(rng.randint(1, 2)): mutate(rng, n, cls)
    rng.shuffle(cls)
    for c in cls: rng.shuffle(c)
    fname = os.path.join(tmpdir, f"color-{seed}.cnf")
    with open(fname, "w") as f:
        f.write(f"p cnf {n} {len(cls)}\n")
        for c in cls: f.write(" ".join(map(str, c)) + " 0\n")
    opts = ["--colorcnt", "1", "--colornodes", str(rng.choice([0, 10, 1000, 250000])),
            "--colordomnodes", str(rng.choice([0, 0, 5, 1000, 100000]))]
    got, info, out = run(ganak, opts, fname, timeout)
    ref, _, _ = run(ganak, ["--colorcnt", "0"], fname, timeout)
    msg = None
    if got == "NOCOUNT" or ref == "NOCOUNT":
        msg = f"no count / crash: {out[-500:]}"
    elif got is not None and ref is not None and got != ref:
        msg = f"WRONG: color {got} vs normal {ref}"
    elif not mutated and info is None and got is not None and 0 < len(cls) <= 64000:
        msg = "unmutated coloring CNF was not recognized"
    if msg is None:
        os.unlink(fname)
        return seed, None, info is not None and "counted" in info.group(0), info is not None
    return seed, f"{msg}\n  ganak {' '.join(opts)} {fname}", False, info is not None

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--num", type=int, default=300)
    p.add_argument("--ganak", default="build/ganak")
    p.add_argument("--seed", type=int, default=random.randrange(1 << 30))
    p.add_argument("-j", type=int, default=16)
    p.add_argument("-t", type=float, default=20)
    a = p.parse_args()
    tmpdir = "fuzz-color-out"
    os.makedirs(tmpdir, exist_ok=True)
    fails = counted = decoded = 0
    with Pool(a.j) as pool:
        for seed, err, cnt, dec in pool.imap_unordered(one, [(a.seed+i, a.ganak, tmpdir, a.t) for i in range(a.num)]):
            counted += cnt
            decoded += dec
            if err:
                fails += 1
                print(f"seed {seed}: {err}", flush=True)
    print(f"{a.num} runs, {decoded} recognized, {counted} counted directly, {fails} failures (base seed {a.seed})")
    sys.exit(1 if fails else 0)

if __name__ == "__main__":
    main()
