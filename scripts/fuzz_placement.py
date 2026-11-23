#!/usr/bin/env python3
# Fuzz the placement counter (--placementcnt) against normal counting.
# Random Sudoku-like labeled exact covers (Latin squares, Sudokus with any box
# shape, diagonal and irregular-region variants, cells that exclude some
# labels), in the pairwise and the capacity encoding, with random clues,
# renaming, and mutations that must make the counter decline or still agree.
# Usage: scripts/fuzz_placement.py [--num N] [--ganak build/ganak] [-j 16]
import argparse, os, random, re, subprocess, sys
from multiprocessing.pool import ThreadPool as Pool  # workers only wait on subprocesses

def regions(rng, n):
    """Constraint families over the n*n grid, each a list of cell lists."""
    rows = [[r*n + c for c in range(n)] for r in range(n)]
    cols = [[r*n + c for r in range(n)] for c in range(n)]
    fams = [rows, cols]
    shapes = [(a, n // a) for a in range(2, n) if n % a == 0]
    kind = rng.randrange(4)
    if kind >= 1 and shapes:  # boxes
        br, bc = rng.choice(shapes)
        fams.append([[(R*br + i)*n + C*bc + j for i in range(br) for j in range(bc)]
                     for R in range(n // br) for C in range(n // bc)])
    if kind == 3:  # diagonals
        fams.append([[i*n + i for i in range(n)]])
        fams.append([[i*n + n-1-i for i in range(n)]])
    return fams

def solve_one(rng, n, fams, allowed, limit=20000):
    """Random solution by backtracking, or None (keeps the clues consistent)."""
    grid = [-1] * (n*n)
    member = [[] for _ in range(n*n)]
    for fam in fams:
        for g in fam:
            for c in g: member[c].append(tuple(g))
    used = {}
    steps = [0]
    def ok(c, d):
        return allowed[c][d] and all(d not in used.get(g, ()) for g in member[c])
    def go(i):
        steps[0] += 1
        if steps[0] > limit: return False
        if i == n*n: return True
        ds = list(range(n)); rng.shuffle(ds)
        for d in ds:
            if ok(i, d):
                grid[i] = d
                for g in member[i]: used.setdefault(g, set()).add(d)
                if go(i+1): return True
                for g in member[i]: used[g].discard(d)
        grid[i] = -1
        return False
    return grid if go(0) else None

def build(rng):
    n = rng.choice([3, 4, 4, 4, 5, 6, 6, 9])
    fams = regions(rng, n)
    allowed = [[True]*n for _ in range(n*n)]
    if rng.random() < 0.2:  # some cells exclude a label (no variable for it)
        for _ in range(rng.randint(1, n)):
            allowed[rng.randrange(n*n)][rng.randrange(n)] = False
    sol = solve_one(rng, n, fams, allowed)
    # Clue density: 9x9 needs many clues for normal counting to finish
    lo = 0.55 if n == 9 else (0.0 if n <= 5 else 0.25)
    clues = []
    if sol is not None:
        p = rng.uniform(lo, min(1.0, lo + 0.4))
        clues = [(c, sol[c]) for c in range(n*n) if rng.random() < p]
    elif n == 9:
        return None
    if rng.random() < 0.15 and clues:  # a clue that may contradict
        c = rng.randrange(n*n); clues.append((c, rng.randrange(n)))
    neg = []
    if rng.random() < 0.2:
        neg = [(rng.randrange(n*n), rng.randrange(n)) for _ in range(rng.randint(1, 3))]
    # variables
    var = {}
    for c in range(n*n):
        for d in range(n):
            if allowed[c][d]: var[c, d] = len(var) + 1
    scopes = [[var[c, d] for d in range(n) if (c, d) in var] for c in range(n*n)]  # cells
    other = []
    for fam in fams:
        for g in fam:
            for d in range(n):
                s = [var[c, d] for c in g if (c, d) in var]
                other.append(s)
    if any(len(s) == 0 for s in scopes + other): return None
    # A one-variable constraint is written as a unit clause, so it cannot serve
    # as a cell or label constraint and the CNF may legitimately be declined
    degenerate = any(len(s) == 1 for s in scopes + other)
    cls = []
    # The capacity certificate needs every variable in equally many constraints
    degree = {c: sum(c in g for fam in fams for g in fam) for c in range(n*n)}
    capacity = (rng.random() < 0.4 and len(set(degree.values())) == 1
                and all(all(row) for row in allowed))
    if capacity:
        # cells as cliques, other constraints at-least-one only
        for s in scopes:
            cls.append(list(s))
            cls += [[-a, -b] for i, a in enumerate(s) for b in s[i+1:]]
        cls += [list(s) for s in other]
    else:
        for s in scopes + other:
            cls.append(list(s))
            cls += [[-a, -b] for i, a in enumerate(s) for b in s[i+1:]]
    for c, d in clues:
        if (c, d) in var: cls.append([var[c, d]])
    for c, d in neg:
        if (c, d) in var: cls.append([-var[c, d]])
    return len(var), cls, n, degenerate

def mutate(rng, nv, cls):
    r = rng.randrange(6)
    if r == 0: del cls[rng.randrange(len(cls))]
    elif r == 1: cls.append([rng.choice([-1, 1]) * rng.randint(1, nv)])
    elif r == 2:
        c = cls[rng.randrange(len(cls))]
        c[rng.randrange(len(c))] *= -1
    elif r == 3:
        a, b = rng.sample(range(1, nv+1), 2)
        cls.append([-a, -b])
    elif r == 4:
        c = cls[rng.randrange(len(cls))]
        if len(c) > 1: c.pop(rng.randrange(len(c)))
    else:
        cls.append(list(cls[rng.randrange(len(cls))]))

def run(ganak, opts, fname, timeout):
    try:
        out = subprocess.run([ganak] + opts + [fname], capture_output=True, text=True, timeout=timeout).stdout
    except subprocess.TimeoutExpired:
        return None, None, ""
    m = re.search(r"^c s exact arb int (\d+)", out, re.M)
    return (int(m.group(1)) if m else "NOCOUNT"), re.search(r"\[placement\].*", out), out

def one(args):
    seed, ganak, tmpdir, timeout = args
    rng = random.Random(seed)
    built = None
    while built is None: built = build(rng)
    nv, cls, n, degenerate = built
    mutated = rng.random() < 0.3
    if mutated:
        for _ in range(rng.randint(1, 2)): mutate(rng, nv, cls)
    perm = list(range(1, nv+1)); rng.shuffle(perm)
    cls = [[(1 if l > 0 else -1) * perm[abs(l)-1] for l in c] for c in cls]
    rng.shuffle(cls)
    for c in cls: rng.shuffle(c)
    fname = os.path.join(tmpdir, f"placement-{seed}.cnf")
    with open(fname, "w") as f:
        f.write(f"p cnf {nv} {len(cls)}\n")
        for c in cls: f.write(" ".join(map(str, c)) + " 0\n")
    opts = ["--placementcnt", "1",
            "--placementbudget", str(rng.choice([0, 3, 100, 10000, 2000000000, 2000000000])),
            "--placementdomain", str(rng.choice([1, 5, 200000, 200000]))]
    got, info, out = run(ganak, opts, fname, timeout)
    ref, _, _ = run(ganak, ["--placementcnt", "0"], fname, timeout)
    msg = None
    if got == "NOCOUNT" or ref == "NOCOUNT":
        msg = f"no count / crash: {out[-500:]}"
    elif got is not None and ref is not None and got != ref:
        msg = f"WRONG: placement {got} vs normal {ref}"
    elif not mutated and not degenerate and info is None and got is not None:
        msg = "unmutated labeled exact cover was not recognized"
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
    p.add_argument("-t", type=float, default=30)
    a = p.parse_args()
    tmpdir = "fuzz-placement-out"
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
