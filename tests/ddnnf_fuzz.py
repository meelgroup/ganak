#!/usr/bin/env python3
"""Fuzz-test Ganak's d-DNNF compilation (`ganak --compile`).

For each random CNF:
  * brute-force the exact model count and model set (oracle),
  * compile the CNF to a d4 .nnf with `ganak --compile`,
  * check the circuit's structural model count == brute count (STRONG mode),
  * check the circuit's model set == brute model set (validates arc literals),
  * also sanity-check against Ganak's own reported count.

For --weak it only checks the file is well-formed and counts without crashing
(the weak count is intentionally wrong), and reports how much smaller the weak
circuit is.

Usage: ddnnf_fuzz.py [num_tests] [--weak] [--seed N]
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


def compile_nnf(path, nnf, weak_level):
    args = [GANAK, "--compile", nnf]
    if weak_level:
        args += ["--weak", str(weak_level)]
    args.append(path)
    r = subprocess.run(args, capture_output=True, text=True)
    return r


def parse_monotone(stdout):
    for line in stdout.splitlines():
        if line.startswith("c o [compile-weak] monotone vars:"):
            return set(int(x) for x in line.split(":", 1)[1].split())
    return set()


def main():
    n = 200
    weak_level = 0   # 0=strong, 1=global monotone, 2=residual monotone
    seed = random.randrange(1 << 30)
    args = sys.argv[1:]
    i = 0
    while i < len(args):
        if args[i] == "--weak":
            if i + 1 < len(args) and args[i + 1].isdigit():
                weak_level = int(args[i + 1])
                i += 1
            else:
                weak_level = 1
        elif args[i] == "--seed":
            seed = int(args[i + 1])
            i += 1
        else:
            n = int(args[i])
        i += 1
    weak = weak_level > 0
    random.seed(seed)
    print(f"seed={seed} tests={n} mode={'WEAK-' + str(weak_level) if weak else 'STRONG'}")

    fails = 0
    weak_smaller = 0
    weak_total = 0
    weak_not_wdnnf = 0
    weak_faithful = 0

    # Reserve race-proof temp paths for this process (atomic O_CREAT|O_EXCL).
    cnf = unique_file("fz", ".cnf")
    nnf = cnf + ".nnf"
    strong_nnf = cnf + ".strong.nnf"

    def save_fail(t):
        shutil.copy(cnf, unique_file("fail", ".cnf"))
        if os.path.exists(nnf):
            shutil.copy(nnf, unique_file("fail", ".nnf"))

    for t in range(n):
        nv = random.randint(7, 16)
        nc = random.randint(nv, nv * 4)
        k = random.choice([2, 3, 3, 4])
        clauses = gen_cnf(nv, nc, k)
        write_cnf(cnf, nv, clauses)
        bc, bmodels = brute(nv, clauses)

        if os.path.exists(nnf):
            os.remove(nnf)
        r = compile_nnf(cnf, nnf, weak_level)
        if not os.path.exists(nnf):
            print(f"FAIL[{t}] no .nnf produced. stderr:\n{r.stderr}")
            save_fail(t)
            fails += 1
            continue
        nodes, arcs, root = dv.parse(nnf)
        sc = dv.count(nodes, arcs, root)

        # strong circuit node count for comparison
        compile_nnf(cnf, strong_nnf, 0)
        snodes, sarcs, sroot = dv.parse(strong_nnf)

        if weak:
            weak_total += 1
            if len(nodes) < len(snodes):
                weak_smaller += 1
            if sc < 0:
                print(f"FAIL[{t}] weak count negative {sc}")
                fails += 1
                continue
            # FAITHFULNESS (what matters for functional synthesis): evaluated as
            # a Boolean FUNCTION on complete assignments, does the circuit equal F?
            # The current weak cut is an approximation (it drops/breaks clause
            # constraints), so this is INFORMATIONAL only -- weak is wrong by design.
            fmodels = dv.function_models(nodes, arcs, root, nv)
            if fmodels == bmodels:
                weak_faithful += 1
            ok_w, _ = dv.check_weak_decomposable(nodes, arcs, root)
            if not ok_w:
                weak_not_wdnnf += 1
            continue

        ok, msg = dv.check_decomposable(nodes, arcs, root)
        if not ok:
            print(f"FAIL[{t}] strong circuit not decomposable: {msg}")
            save_fail(t)
            fails += 1
            continue

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
        # FAITHFUL AS A FUNCTION: evaluated on every complete assignment, the
        # circuit must equal F. This is exactly the property functional synthesis
        # needs, so we assert it for the strong (d-DNNF) circuit.
        fmodels = dv.function_models(nodes, arcs, root, nv)
        if fmodels != bmodels:
            print(f"FAIL[{t}] strong circuit not faithful as a function "
                  f"|circuit|={len(fmodels)} |brute|={len(bmodels)}")
            save_fail(t)
            fails += 1
            continue

    if weak:
        print(f"done {n} weak tests, {fails} failures; "
              f"FAITHFUL-as-function (synthesis-usable): {weak_faithful}/{weak_total}; "
              f"NOT weak-decomposable: {weak_not_wdnnf}/{weak_total}; "
              f"smaller than strong: {weak_smaller}/{weak_total}")
    else:
        print(f"done {n} strong tests, {fails} failures")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
