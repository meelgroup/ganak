#!/usr/bin/env python3
"""Functional-synthesis round-trip test for `ganak --compile`.

For each random projected CNF (inputs X = sampling vars 1..k, outputs Y = rest),
compile to a circuit; for every satisfiable X, synthesize psi(X) off the circuit
and check F(X, psi(X)) holds. The real correctness criterion for synthesis (not
model counting or faithfulness).

Run `ddnnf_synth.py --help` for all options.
"""
import argparse
import itertools
import os
import random
import shutil
import stat
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ddnnf_verify as dv

GANAK = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "build", "ganak")
TMP = "/tmp/ddnnf_synth"
os.makedirs(TMP, exist_ok=True)


def unique_file(begin, end=".cnf"):
    i = 1
    while True:
        f = os.path.join(TMP, f"{begin}_{i}{end}")
        try:
            fd = os.open(f, os.O_CREAT | os.O_EXCL, stat.S_IREAD | stat.S_IWRITE)
            os.close(fd)
            return f
        except OSError:
            i += 1


def gen(nv, nc):
    cls = []
    for _ in range(nc):
        c = set()
        while len(c) < min(3, nv):
            v = random.randint(1, nv)
            c.add(v if random.random() < 0.5 else -v)
        cls.append(list(c))
    return cls


def write_cnf(path, nv, k, cls):
    with open(path, "w") as f:
        f.write("c p show " + " ".join(str(v) for v in range(1, k + 1)) + " 0\n")
        f.write(f"p cnf {nv} {len(cls)}\n")
        for c in cls:
            f.write(" ".join(map(str, c)) + " 0\n")


def sat(cls, a):
    return all(any((l > 0) == a[abs(l)] for l in c) for c in cls)


def parse_args(argv):
    p = argparse.ArgumentParser(
        prog="ddnnf_synth.py",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--num", type=int, default=None, metavar="N",
                   help="number of random instances to test (default: run forever)")
    p.add_argument("--synthesis", action="store_true",
                   help="compile with --synthesis 1 (wDNNF share-and-branch) and "
                        "report circuit size vs the faithful compile")
    p.add_argument("--seed", type=int, default=None,
                   help="RNG seed (default: random; printed so runs can be reproduced)")
    p.add_argument("--satoff", action="store_true",
                   help="pass --satsolver 0 to ganak (disable the SAT oracle)")
    p.add_argument("--minvars", type=int, default=6, metavar="N",
                   help="min variables per instance (default: 6)")
    p.add_argument("--maxvars", type=int, default=14, metavar="N",
                   help="max variables per instance; the brute-force oracle is 2^N, "
                        "so this is the main runtime/coverage knob (default: 14)")
    p.add_argument("-v", "--verbose", action="store_true",
                   help="print per-instance progress (params, #solvable X, OK/FAIL)")
    p.add_argument("-q", "--quiet", action="store_true",
                   help="only print the final summary line")
    p.add_argument("--keep-going", action="store_true",
                   help="run all instances instead of stopping at the first wrong "
                        "case, and tally every failing X within each instance "
                        "(default: stop immediately on the first failure)")
    args = p.parse_args(argv)
    if args.minvars < 1 or args.maxvars < args.minvars:
        p.error("need 1 <= --minvars <= --maxvars")
    if args.seed is None:
        args.seed = random.randrange(1 << 30)
    return args


def main():
    args = parse_args(sys.argv[1:])
    n = args.num   # None => run forever (until a failure, with fail-fast)
    synth = args.synthesis
    random.seed(args.seed)
    if not args.quiet:
        print(f"seed={args.seed} tests={n if n is not None else 'forever'} "
              f"synthesis={synth} vars=[{args.minvars}..{args.maxvars}]")

    cnf = unique_file("s", ".cnf")
    nnf = cnf + ".nnf"
    fails = 0
    no_witness = 0
    size_w = 0
    size_0 = 0
    ran = 0
    t = 0
    while n is None or t < n:
        t += 1
        ran += 1
        nv = random.randint(args.minvars, args.maxvars)
        k = random.randint(2, max(2, nv - 2))   # |X|
        nc = random.randint(nv, nv * 3)
        cls = gen(nv, nc)
        write_cnf(cnf, nv, k, cls)

        if os.path.exists(nnf):
            os.remove(nnf)
        a = [GANAK, "--compile", nnf]
        if synth:
            a += ["--synthesis", "1"]
        if args.satoff:
            a += ["--satsolver", "0"]
        a.append(cnf)
        r = subprocess.run(a, capture_output=True, text=True)
        if not os.path.exists(nnf):
            print(f"FAIL[{t}] no .nnf. stderr:\n{r.stderr[-500:]}")
            shutil.copy(cnf, unique_file("fail", ".cnf"))
            fails += 1
            if not args.keep_going:
                break
            continue
        nodes, arcs, root = dv.parse(nnf)
        size_w += len(nodes)
        if synth:  # compile faithful (no --synthesis) too, to compare circuit size
            n0 = nnf + ".w0"
            subprocess.run([GANAK, "--compile", n0, cnf], capture_output=True)
            if os.path.exists(n0):
                nn0, _, _ = dv.parse(n0)
                size_0 += len(nn0)

        # which input assignments X are satisfiable?
        solvable = set()
        for bits in itertools.product((False, True), repeat=nv):
            assign = {v + 1: bits[v] for v in range(nv)}
            if sat(cls, assign):
                solvable.add(tuple(bits[:k]))

        bad = 0
        for xb in solvable:
            xassign = {v + 1: xb[v] for v in range(k)}
            psi = dv.synthesize(nodes, arcs, root, xassign)
            if psi is None:
                no_witness += 1
                bad += 1
                if not args.keep_going:
                    break
                continue
            full = {v: psi.get(v, False) for v in range(1, nv + 1)}
            full.update(xassign)
            if not sat(cls, full):
                bad += 1
                if not args.keep_going:
                    break
        if bad:
            print(f"FAIL[{t}] synthesis invalid (nv={nv} k={k} nc={nc} "
                  f"bad_X={bad}/{len(solvable)})")
            shutil.copy(cnf, unique_file("fail", ".cnf"))
            shutil.copy(nnf, unique_file("fail", ".nnf"))
            fails += 1
            if not args.keep_going:
                print(f"stopping at first failure (seed={args.seed}); "
                      f"saved to {TMP}/fail_*.{{cnf,nnf}}. Use --keep-going to run all.")
                break
        elif args.verbose:
            print(f"ok[{t}] nv={nv} k={k} nc={nc} nodes={len(nodes)} "
                  f"solvable_X={len(solvable)}")

    msg = (f"done {ran} synth tests, {fails} failures "
           f"({no_witness} satisfiable-X with no witness)")
    if synth and size_0:
        msg += f"; circuit size synthesis/faithful = {size_w/size_0:.2f}"
    print(msg)
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
