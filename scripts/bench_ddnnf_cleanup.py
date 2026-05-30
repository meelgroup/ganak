#!/usr/bin/env python3
"""Performance + regression bench for `build/ddnnf-cleanup`.

What this does
--------------
1. Generates random uniform 3-CNFs at sizes nv in [40, 70, 80, 90, 100] (configurable),
   nc = 3*nv. CNFs go to /tmp and are NOT committed -- they regenerate each run.
2. Compiles each with `ganak --compile` to a .nnf.
3. Runs `ddnnf-cleanup` twice per CNF: default (= strict-decomp ON) and
   `--no-strict-decomp` (= just BFS-renumber + drop dead nodes).
4. Records best-of-3 wall time and peak RSS for each, plus output file size,
   number of strict-decomp fixes applied (parsed from stderr), and a
   correctness cross-check against ganak's standalone count.
5. Also runs the seed=207900131-derived regression CNF (the smoking-gun
   instance where the raw circuit isn't strict-decomposable and the cleanup's
   strict-decomp pass repairs it).

When to run
-----------
NOT part of the ctest suite (would dominate wall time at large sizes). Run it
by hand when:
  - touching the strict-decomp algorithm (compute_subtree_vars, compute_forced,
    scrub_var, cleanup_decomp) in src/ddnnf_cleanup.cpp;
  - changing the bitset representation (VarSet/LitSet);
  - investigating a memory blow-up on a large .nnf.

Usage:  ./scripts/bench_ddnnf_cleanup.py [--sizes 50,70,80,90,100] [--reps 3]
                                         [--timeout 180] [--seed 1]
"""
import argparse
import os
import random
import re
import shutil
import subprocess
import sys
import time


ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GANAK = os.path.join(ROOT, "build", "ganak")
CLEANUP = os.path.join(ROOT, "build", "ddnnf-cleanup")
VERIFY = os.path.join(ROOT, "tests", "ddnnf_verify.py")
REGRESSION_CNF = os.path.join(ROOT, "tests", "cnf-files", "ddnnf",
                              "known_nondecomp_regression.cnf")
TMP = "/tmp/bench_ddnnf_cleanup"


def gen_cnf(path, nv, nc, k, rng):
    with open(path, "w") as f:
        f.write(f"p cnf {nv} {nc}\n")
        for _ in range(nc):
            cl = set()
            while len(cl) < min(k, nv):
                v = rng.randint(1, nv)
                cl.add(v if rng.random() < 0.5 else -v)
            f.write(" ".join(map(str, cl)) + " 0\n")


def time_proc(argv, timeout):
    """Run argv under /usr/bin/time, return (wall_seconds, peak_kb, stderr, ok)."""
    full = ["/usr/bin/time", "-f", "%e %M"] + argv
    try:
        r = subprocess.run(full, capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        return (timeout, 0, "TIMEOUT", False)
    # /usr/bin/time prints "<wall> <kb>" on its own line at the end of stderr
    m = re.search(r"(\d+(?:\.\d+)?)\s+(\d+)\s*$", r.stderr.strip())
    if not m:
        return (None, None, r.stderr, False)
    return (float(m.group(1)), int(m.group(2)), r.stderr, r.returncode == 0)


def best_of(argv, reps, timeout):
    """Take min wall + corresponding mem across reps."""
    best = (float("inf"), 0, "", False)
    for _ in range(reps):
        t, m, err, ok = time_proc(argv, timeout)
        if t is None:
            return (None, 0, err, False)
        if t < best[0]:
            best = (t, m, err, ok)
    return best


def ganak_count(cnf):
    r = subprocess.run([GANAK, cnf], capture_output=True, text=True, timeout=300)
    for line in r.stdout.splitlines():
        if "exact arb" in line:
            return line.strip().split()[-1]
    return None


def verify_count(nnf):
    r = subprocess.run(["python3", VERIFY, nnf], capture_output=True, text=True, timeout=300)
    return r.stdout.strip().splitlines()[-1] if r.stdout.strip() else None


def parse_fixes(stderr):
    m = re.search(r"strict-decomp iters=(\d+) fixes=(\d+)", stderr or "")
    return (int(m.group(1)), int(m.group(2))) if m else (0, 0)


def fmt_mb(kb):
    return "—" if not kb else f"{kb/1024:.0f}"


def fmt_t(t):
    if t is None:
        return "FAIL"
    if t > 60:
        return f"{t:.0f}s"
    return f"{t:.2f}s"


def run_one(name, cnf, nnf, args):
    """Compile + bench both modes. Returns a dict row."""
    if os.path.exists(nnf):
        os.remove(nnf)
    subprocess.run([GANAK, "--compile", nnf, cnf], capture_output=True, timeout=300)
    if not os.path.exists(nnf):
        print(f"  {name}: ganak compile FAILED", file=sys.stderr)
        return None
    nodes = sum(1 for line in open(nnf) if line[:1] in "oatf")
    in_size = os.path.getsize(nnf)

    # warmup so file is in page cache
    with open(nnf, "rb") as f: f.read()
    # Per-case scratch paths so a previous case's leftover file never gets
    # mistaken for this case's output if the cleanup gets OOM-killed.
    out_bfs = nnf + ".bfs.out"
    out_str = nnf + ".str.out"
    for p in (out_bfs, out_str):
        if os.path.exists(p): os.remove(p)
    subprocess.run([CLEANUP, "--no-strict-decomp", nnf, out_bfs], capture_output=True)
    subprocess.run([CLEANUP, nnf, out_str], capture_output=True, timeout=args.timeout)

    # Re-delete before timing so each rep starts with a known state. If a rep
    # crashes, best_of will see no output file and report it via the count
    # mismatch path below.
    for _ in range(args.reps):
        pass  # no-op
    bfs = best_of([CLEANUP, "--no-strict-decomp", nnf, out_bfs], args.reps, args.timeout)
    strc = best_of([CLEANUP, nnf, out_str], args.reps, args.timeout)
    iters, fixes = parse_fixes(strc[2])

    expected = ganak_count(cnf)
    # Correctness: cleaned-strict file's count should match ganak. Existence
    # check guards against an OOM-killed cleanup leaving us with no output.
    if not os.path.exists(out_str) or os.path.getsize(out_str) == 0:
        match = f"FAIL(no output, ganak={expected})"
    else:
        str_cnt = verify_count(out_str)
        match = "OK" if str_cnt == expected else f"MISMATCH({str_cnt} vs {expected})"

    return {
        "name": name,
        "nodes": nodes,
        "in_mb": in_size / 1024 / 1024,
        "bfs_t": bfs[0], "bfs_mb": bfs[1],
        "str_t": strc[0], "str_mb": strc[1],
        "fixes": fixes, "iters": iters,
        "match": match,
    }


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--sizes", default="40,70,80,90,100",
                   help="comma-separated nv values; nc = nv*3 (default: %(default)s)")
    p.add_argument("--reps", type=int, default=3, help="best-of-N timing (default: 3)")
    p.add_argument("--timeout", type=int, default=180,
                   help="per-cleanup-invocation timeout in seconds (default: 180)")
    p.add_argument("--seed", type=int, default=1, help="RNG seed for CNF generation")
    p.add_argument("--keep", action="store_true",
                   help="don't delete the generated /tmp CNFs + .nnf outputs on exit")
    args = p.parse_args()

    if not (os.path.exists(GANAK) and os.path.exists(CLEANUP)):
        print(f"ERROR: build/ganak or build/ddnnf-cleanup not found. Build first.",
              file=sys.stderr)
        sys.exit(1)

    os.makedirs(TMP, exist_ok=True)
    rng = random.Random(args.seed)
    sizes = [int(s) for s in args.sizes.split(",") if s.strip()]

    rows = []

    # Regression instance first (small, fast, repairs a real violation).
    reg_cnf = os.path.join(TMP, "regression.cnf")
    shutil.copy(REGRESSION_CNF, reg_cnf)
    nnf = os.path.join(TMP, "regression.nnf")
    row = run_one("regression(16v/40c)", reg_cnf, nnf, args)
    if row:
        rows.append(row)

    # Random 3-CNFs at the configured sizes.
    for nv in sizes:
        nc = nv * 3
        cnf = os.path.join(TMP, f"random_nv{nv}.cnf")
        gen_cnf(cnf, nv, nc, 3, rng)
        nnf = os.path.join(TMP, f"random_nv{nv}.nnf")
        row = run_one(f"nv={nv}/nc={nc}", cnf, nnf, args)
        if row:
            rows.append(row)

    # Print a table.
    headers = ["instance", "nodes", "in MB", "bfs t", "bfs MB",
               "strict t", "strict MB", "fixes/iters", "count"]
    widths = [22, 10, 8, 9, 8, 9, 9, 12, 10]
    def fmt_row(cells):
        return "  ".join(f"{c:<{w}}" for c, w in zip(cells, widths))

    print()
    print(fmt_row(headers))
    print(fmt_row(["-" * (w - 1) for w in widths]))
    for r in rows:
        print(fmt_row([
            r["name"],
            f"{r['nodes']:,}",
            f"{r['in_mb']:.1f}",
            fmt_t(r["bfs_t"]),
            fmt_mb(r["bfs_mb"]),
            fmt_t(r["str_t"]),
            fmt_mb(r["str_mb"]),
            f"{r['fixes']}/{r['iters']}",
            r["match"],
        ]))
    print()
    print(f"seed={args.seed} reps={args.reps} timeout={args.timeout}s")
    print(f"(temp files in {TMP}, --keep to retain)")

    if not args.keep:
        shutil.rmtree(TMP, ignore_errors=True)

    if any(r["match"] != "OK" for r in rows):
        sys.exit(2)


if __name__ == "__main__":
    main()
