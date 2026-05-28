#!/usr/bin/env python3
"""Functional-synthesis round-trip test for `ganak --compile`.

For each random projected CNF (inputs X = sampling vars 1..k, outputs Y = the
rest), compile to a circuit, then for every input assignment X that the formula
can satisfy, synthesize an output assignment psi(X) by reading it off the
circuit and check that F(X, psi(X)) actually holds. This is the real
correctness criterion for using the (weak) circuit in Boolean functional
synthesis -- it does not care about model counting or faithfulness-as-a-function.

Usage: ddnnf_synth.py [num_tests] [--weak N] [--seed S]
"""
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


def gen(nv, nc, k):
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


def main():
    n = 200
    weak = 0
    satoff = False
    seed = random.randrange(1 << 30)
    args = sys.argv[1:]
    i = 0
    while i < len(args):
        if args[i] == "--weak":
            weak = int(args[i + 1]); i += 1
        elif args[i] == "--seed":
            seed = int(args[i + 1]); i += 1
        elif args[i] == "--satoff":
            satoff = True
        else:
            n = int(args[i])
        i += 1
    random.seed(seed)
    print(f"seed={seed} tests={n} weak={weak}")

    cnf = unique_file("s", ".cnf")
    nnf = cnf + ".nnf"
    fails = 0
    no_witness = 0
    for t in range(n):
        nv = random.randint(6, 14)
        k = random.randint(2, max(2, nv - 2))   # |X|
        nc = random.randint(nv, nv * 3)
        cls = gen(nv, nc, 3)
        write_cnf(cnf, nv, k, cls)

        if os.path.exists(nnf):
            os.remove(nnf)
        a = [GANAK, "--compile", nnf]
        if weak:
            a += ["--weak", str(weak)]
        if satoff:
            a += ["--satsolver", "0"]
        a.append(cnf)
        r = subprocess.run(a, capture_output=True, text=True)
        if not os.path.exists(nnf):
            print(f"FAIL[{t}] no .nnf. stderr:\n{r.stderr[-500:]}")
            fails += 1
            continue
        nodes, arcs, root = dv.parse(nnf)

        # which input assignments X are satisfiable (extendable to a full model)?
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
                break
            full = {v: psi.get(v, False) for v in range(1, nv + 1)}
            full.update(xassign)
            if not sat(cls, full):
                bad += 1
                break
        if bad:
            print(f"FAIL[{t}] synthesis invalid (nv={nv} k={k} nc={nc})")
            shutil.copy(cnf, unique_file("fail", ".cnf"))
            shutil.copy(nnf, unique_file("fail", ".nnf"))
            fails += 1

    print(f"done {n} synth tests, {fails} failures "
          f"({no_witness} with a satisfiable X that the circuit gave no witness for)")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
