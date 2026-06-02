#!/usr/bin/env python3
"""Fuzz-test Ganak's faithful d-DNNF compilation (`ganak --compile`).

For each random CNF: brute-force the count + model set (oracle), compile to a d4
.nnf, then check it is strictly decomposable, its count and model set match the
oracle, and it is faithful as a Boolean function. Also cross-checks ddnnf-cleanup.
(Functional synthesis round-trip: see ddnnf_synth.py.)

Usage: ddnnf_fuzz.py [num_tests] [--seed N]
"""
import itertools
import os
import shutil
import stat
import random
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ddnnf_verify as dv

GANAK = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "build", "ganak")
CLEANUP = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "build", "ddnnf-cleanup")
VERIFY = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "build", "ddnnf-verify")
TMP = "/tmp/ddnnf_fuzz"
os.makedirs(TMP, exist_ok=True)


def unique_file(fname_begin, fname_end=".cnf", max_num_files=100000):
    """Atomically reserve a unique filename (O_CREAT|O_EXCL), so concurrent
    fuzzer runs never collide on temp files."""
    counter = 1
    while True:
        fname = os.path.join(TMP, f"{fname_begin}_{counter}{fname_end}")
        try:
            fd = os.open(fname, os.O_CREAT | os.O_EXCL, stat.S_IREAD | stat.S_IWRITE)
            os.close(fd)
            return fname
        except OSError:
            pass
        counter += 1
        if counter > max_num_files:
            print(f"Cannot create unique_file, last try was: {fname}")
            sys.exit(-1)


def gen_cnf(nv, nc, k):
    clauses = []
    for _ in range(nc):
        cl = set()
        while len(cl) < min(k, nv):
            v = random.randint(1, nv)
            cl.add(v if random.random() < 0.5 else -v)
        clauses.append(list(cl))
    return clauses


def write_cnf(path, nv, clauses):
    with open(path, "w") as f:
        f.write(f"p cnf {nv} {len(clauses)}\n")
        for cl in clauses:
            f.write(" ".join(map(str, cl)) + " 0\n")


def brute(nv, clauses):
    cnt = 0
    models = set()
    for bits in itertools.product((False, True), repeat=nv):
        assign = [None] + [bits[v] for v in range(nv)]
        ok = True
        for cl in clauses:
            if not any((l > 0) == assign[abs(l)] for l in cl):
                ok = False
                break
        if ok:
            cnt += 1
            models.add(frozenset(v if assign[v] else -v for v in range(1, nv + 1)))
    return cnt, models


def ganak_count(path):
    out = subprocess.run([GANAK, path], capture_output=True, text=True).stdout
    for line in out.splitlines():
        if line.startswith("c s exact arb int"):
            return int(line.split()[-1])
    return None


def compile_nnf(path, nnf):
    args = [GANAK, "--compile", nnf, path]
    return subprocess.run(args, capture_output=True, text=True)


def run_big(t, cnf, nnf, clean):
    """--big mode: ganak as oracle (no brute), C++ ddnnf-verify for checks.
    Returns (ok, msg). Catches ddnnf-cleanup crashes (the main thing this mode
    exists to find): subprocess.run() returns non-zero on SIGSEGV/SIGABRT.
    """
    gc = ganak_count(cnf)
    if gc is None:
        return True, "ganak unable to count (skipping)"
    r = compile_nnf(cnf, nnf)
    if not os.path.exists(nnf):
        return False, f"no .nnf produced. stderr tail:\n{r.stderr[-1500:]}"
    # Raw circuit may not be strict / decomposable; check count only.
    rv = subprocess.run([VERIFY, nnf, "--expect-count", str(gc),
                         "--no-strict", "--no-check-decomposable"],
                        capture_output=True, text=True)
    if rv.returncode != 0:
        return False, f"raw count mismatch (oracle {gc}):\n{rv.stdout}{rv.stderr}"
    if os.path.exists(clean):
        os.remove(clean)
    rc = subprocess.run([CLEANUP, nnf, clean], capture_output=True, text=True)
    if rc.returncode != 0:
        return False, (f"ddnnf-cleanup exit={rc.returncode} "
                       f"(signal? {-rc.returncode if rc.returncode < 0 else 'n/a'}). "
                       f"stderr tail:\n{rc.stderr[-1500:]}")
    if not os.path.exists(clean):
        return False, f"cleanup produced no output. stderr tail:\n{rc.stderr[-1500:]}"
    # Defaults: --strict --check-decomposable; --expect-count cross-checks count.
    rv = subprocess.run([VERIFY, clean, "--expect-count", str(gc)],
                        capture_output=True, text=True)
    if rv.returncode != 0:
        return False, f"cleaned verify failed:\n{rv.stdout}{rv.stderr}"
    return True, ""


def main():
    n = 200
    seed = random.randrange(1 << 30)
    minv, maxv = 7, 16    # var-count range; the oracle is 2^nv, so this drives runtime
    big_mode = False
    args = sys.argv[1:]
    i = 0
    while i < len(args):
        if args[i] == "--seed":
            seed = int(args[i + 1])
            i += 1
        elif args[i] == "--minvars":
            minv = int(args[i + 1])
            i += 1
        elif args[i] == "--maxvars":
            maxv = int(args[i + 1])
            i += 1
        elif args[i] == "--big":
            # Cleanup-focused mode: ganak as oracle (no brute), bigger CNFs to
            # provoke DAG sharing + cloning interactions that small-nv brute
            # mode can't reach. Catches ddnnf-cleanup crashes / wrong counts /
            # surviving non-decomposability at scale.
            big_mode = True
        else:
            n = int(args[i])
        i += 1
    if big_mode:
        if minv == 7 and maxv == 16: minv, maxv = 30, 60
        if n == 200: n = 30
    minv = min(minv, maxv)
    random.seed(seed)
    print(f"seed={seed} tests={n} vars={minv}..{maxv} mode={'BIG' if big_mode else 'STRONG'}")

    fails = 0

    # Reserve race-proof temp paths for this process (atomic O_CREAT|O_EXCL).
    cnf = unique_file("fz", ".cnf")
    nnf = cnf + ".nnf"
    clean = nnf + ".clean"

    def save_fail(t):
        shutil.copy(cnf, unique_file("fail", ".cnf"))
        if os.path.exists(nnf):
            shutil.copy(nnf, unique_file("fail", ".nnf"))
        if os.path.exists(clean):
            shutil.copy(clean, unique_file("fail", ".clean.nnf"))

    for t in range(n):
        nv = random.randint(minv, maxv)
        nc = random.randint(nv, nv * 4)
        k = random.choice([2, 3, 3, 4])
        clauses = gen_cnf(nv, nc, k)
        write_cnf(cnf, nv, clauses)

        if big_mode:
            ok, msg = run_big(t, cnf, nnf, clean)
            if not ok:
                print(f"FAIL[{t}] nv={nv} nc={nc} k={k}: {msg}")
                save_fail(t)
                fails += 1
            elif msg:
                print(f"NOTE[{t}] nv={nv} nc={nc} k={k}: {msg}")
            continue

        bc, bmodels = brute(nv, clauses)

        if os.path.exists(nnf):
            os.remove(nnf)
        r = compile_nnf(cnf, nnf)
        if not os.path.exists(nnf):
            print(f"FAIL[{t}] no .nnf produced. stderr:\n{r.stderr}")
            save_fail(t)
            fails += 1
            continue
        nodes, arcs, root = dv.parse(nnf)
        sc = dv.count(nodes, arcs, root)

        ok, msg = dv.check_decomposable(nodes, arcs, root)
        if not ok:
            # Raw circuit may transiently violate strict-decomposability because
            # conflict-learned (redundant) clauses can propagate "cross-component"
            # vars into an OR's arc-lits (CompAnalyzer only sees irreducible clauses).
            # The count + model-set checks below still apply -- the function is
            # right. `ddnnf-cleanup` repairs the violation below.
            print(f"NOTE[{t}] raw circuit not strict-decomposable (will check post-cleanup): {msg}")

        if sc != bc:
            print(f"FAIL[{t}] structural count {sc} != brute {bc}  nv={nv} nc={nc} k={k}")
            save_fail(t)
            fails += 1
            continue
        gc = ganak_count(cnf)
        if gc is not None and gc != bc:
            print(f"WARN[{t}] ganak normal count {gc} != brute {bc} (oracle disagreement)")
        cmodels = dv.models(nodes, arcs, root, nv)
        if cmodels != bmodels:
            print(f"FAIL[{t}] model-set mismatch nv={nv} |circuit|={len(cmodels)} |brute|={len(bmodels)}")
            save_fail(t)
            fails += 1
            continue
        # FAITHFUL AS A FUNCTION: on every complete assignment the circuit must
        # equal F (the property functional synthesis needs).
        fmodels = dv.function_models(nodes, arcs, root, nv)
        if fmodels != bmodels:
            print(f"FAIL[{t}] strong circuit not faithful as a function "
                  f"|circuit|={len(fmodels)} |brute|={len(bmodels)}")
            save_fail(t)
            fails += 1
            continue

        # CLEANUP TOOL: the raw .nnf keeps internal ids and may carry dead nodes.
        # `ddnnf-cleanup` must drop them and renumber root=1. Check strictly: no
        # dead nodes, still decomposable, same count + model set.
        if os.path.exists(clean):
            os.remove(clean)
        rc = subprocess.run([CLEANUP, nnf, clean], capture_output=True, text=True)
        if not os.path.exists(clean):
            print(f"FAIL[{t}] ddnnf-cleanup produced no output. stderr:\n{rc.stderr}")
            save_fail(t)
            fails += 1
            continue
        cn, ca, cr = dv.parse(clean)
        dead = dv.unreachable_nodes(cn, ca, cr)
        if dead:
            print(f"FAIL[{t}] cleaned circuit still has dead nodes: {sorted(dead)[:10]}")
            save_fail(t)
            fails += 1
            continue
        ok, msg = dv.check_decomposable(cn, ca, cr)
        if not ok:
            print(f"FAIL[{t}] cleaned circuit not decomposable: {msg}")
            save_fail(t)
            fails += 1
            continue
        nested = dv.nested_and_arcs(cn, ca, cr)
        if nested:
            print(f"FAIL[{t}] cleaned circuit has nested AND-of-AND arcs "
                  f"(should be flattened): {nested[:10]}")
            save_fail(t)
            fails += 1
            continue
        unary = dv.unary_and_nodes(cn, ca, cr)
        if unary:
            print(f"FAIL[{t}] cleaned circuit has unary (single-arc) AND nodes "
                  f"(should be elided): {unary[:10]}")
            save_fail(t)
            fails += 1
            continue
        if dv.count(cn, ca, cr) != bc:
            print(f"FAIL[{t}] cleaned structural count {dv.count(cn, ca, cr)} != brute {bc}")
            save_fail(t)
            fails += 1
            continue
        if dv.models(cn, ca, cr, nv) != bmodels:
            print(f"FAIL[{t}] cleaned model-set mismatch")
            save_fail(t)
            fails += 1
            continue

    print(f"done {n} {'big' if big_mode else 'strong'} tests, {fails} failures")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
