#!/usr/bin/env python3
"""Checks literal frequencies in sampled models (DIMACS lines) against expected shares.
Usage: sample_marginals.py samples.txt --tol T LIT=P [LIT=P ...]   (P may be a fraction, e.g. 12/13)
"""
import sys
from fractions import Fraction

args = sys.argv[1:]
path = args.pop(0)
tol = 0.02
if args and args[0] == "--tol":
    tol = float(args[1])
    args = args[2:]
models = [set(int(x) for x in l.split()[:-1]) for l in open(path) if l.strip() and not l.startswith("c")]
if not models:
    sys.exit("no samples")
bad = False
for a in args:
    lit, p = a.split("=")
    lit, p = int(lit), float(Fraction(p))
    f = sum(lit in m for m in models) / len(models)
    ok = abs(f - p) <= tol
    bad |= not ok
    print(f"{lit}: {f:.4f} (expected {p:.4f}) {'ok' if ok else 'FAIL'}")
sys.exit(1 if bad else 0)
