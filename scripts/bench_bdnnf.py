#!/usr/bin/env python3
"""Compares bdnnf with ddnnife on competition CNFs shrunk to a few seconds each.
Usage: scripts/bench_bdnnf.py [-j J] [--only NAME ...] [--timeout S]   (~5 min)
"""

import argparse
import gzip
import os
import random
import subprocess
import sys
import tempfile
import time
from concurrent.futures import ProcessPoolExecutor

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SAT = os.path.dirname(ROOT)
GANAK = f"{ROOT}/build/ganak"
BDNNF = f"{ROOT}/build/bdnnf"
DDNNIFE = f"{SAT}/ddnnife/target/release/ddnnife"
CMS = f"{SAT}/cryptominisat/build/cryptominisat5"
MC23 = os.path.expanduser("~/media/counting/mccomp2023/mc2023-track1-mc")

# name, CNF, fraction of the vars fixed to a model's values (shrinks it), queries
# run on it: the full circuits make each query take about a second; 2-wise
# needs shrunk ones (both tools time out on the full ones)
TWISE = ["2-wise"]
REST = ["count", "features", "core", "samples", "count queries", "SAT queries", "atomic sets"]
INSTANCES = [
    ("t015", f"{MC23}/mc2023_track1_015.cnf.gz", 0.0, REST),
    ("t066", f"{MC23}/mc2023_track1_066.cnf.gz", 0.0, REST),
    ("t076", f"{MC23}/mc2023_track1_076.cnf.gz", 0.0, REST),
    ("t084", f"{MC23}/mc2023_track1_084.cnf.gz", 0.0, REST),
    ("t102", f"{MC23}/mc2023_track1_102.cnf.gz", 0.0, REST),
    ("auto1-s", f"{SAT}/ddnnife/example_input/auto1.cnf", 0.01, TWISE),
    ("t015-s", f"{MC23}/mc2023_track1_015.cnf.gz", 0.01, TWISE),
    ("t027-s", f"{MC23}/mc2023_track1_027.cnf.gz", 0.01, TWISE),
    ("t028-s", f"{MC23}/mc2023_track1_028.cnf.gz", 0.03, TWISE),
    ("t066-s", f"{MC23}/mc2023_track1_066.cnf.gz", 0.01, TWISE),
    ("t003-s", f"{MC23}/mc2023_track1_003.cnf.gz", 0.0, TWISE),
]
NQ = 10000  # samples, count queries, SAT queries

# label, bdnnf args, ddnnife args ({c} circuit, {q} query file)
QUERIES = [
    ("count", ["count", "{c}"], ["count"]),
    ("features", ["lits", "{c}"], ["count-features"]),
    ("core", ["core", "{c}"], ["core"]),
    ("samples", ["sample", "{c}", "-n", str(NQ)], ["urs", "-n", str(NQ)]),
    ("count queries", ["batch", "{c}", "{q}"], ["count-queries", "{q}"]),
    ("SAT queries", ["batch", "{c}", "{q}", "--sat"], ["sat", "{q}"]),
    ("atomic sets", ["atomic", "{c}"], ["atomic-sets"]),
    ("2-wise", ["twise", "{c}", "-t", "2"], ["t-wise", "-t", "2"]),
]


def run(cmd, timeout):
    """(seconds, None, stdout), or (None, None, '') on a timeout or failure."""
    t0 = time.time()
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
    try:
        out, _ = p.communicate(timeout=timeout)
    except subprocess.TimeoutExpired:
        p.kill()
        p.communicate()
        return None, None, ""
    if p.returncode != 0:
        return None, None, ""
    return time.time() - t0, None, out


def timed(cmd, timeout):
    """Runs cmd under /usr/bin/time for its peak memory."""
    tf = tempfile.NamedTemporaryFile(delete=False)
    tf.close()
    t, _, out = run(["/usr/bin/time", "-f", "%M", "-o", tf.name] + cmd, timeout)
    mb = None
    if t is not None:
        try:
            mb = int(open(tf.name).read().split()[-1]) / 1024
        except (ValueError, IndexError):
            pass
    os.unlink(tf.name)
    return t, mb, out


def shrink(name, path, frac, work):
    lines = (gzip.open(path, "rt") if path.endswith(".gz") else open(path)).read().split("\n")
    hdr = next(i for i, l in enumerate(lines) if l.startswith("p cnf"))
    nv, nc = int(lines[hdr].split()[2]), int(lines[hdr].split()[3])
    cls = [l for i, l in enumerate(lines) if i != hdr and l.strip() and not l.startswith(("c", "p"))]
    units = []
    if frac > 0:
        open(f"{work}/{name}.orig.cnf", "w").write("\n".join([lines[hdr]] + cls) + "\n")
        out = subprocess.run([CMS, "--verb", "0", "--polar", "rnd", "--seed", "1", f"{work}/{name}.orig.cnf"],
                             capture_output=True, text=True).stdout
        model = [int(x) for l in out.split("\n") if l.startswith("v") for x in l.split()[1:] if x != "0"]
        units = random.Random(1).sample(model, int(frac * len(model)))
    cnf = f"{work}/{name}.cnf"
    open(cnf, "w").write(f"p cnf {nv} {nc + len(units)}\n" + "\n".join(cls + [f"{u} 0" for u in units]) + "\n")
    return cnf, nv


def one(job):
    name, path, frac, labels, work, tl = job
    cnf, nv = shrink(name, path, frac, work)
    bd, nnf, q = f"{work}/{name}.bdnnf", f"{work}/{name}.nnf", f"{work}/{name}.q"
    tc, mc, _ = timed([GANAK, "--verb", "0", "--compile", bd, cnf], 600)
    subprocess.run([BDNNF, "to-d4", bd, nnf], capture_output=True)
    rng = random.Random(2)
    with open(q, "w") as f:
        for _ in range(NQ):
            f.write(" ".join(str(v if rng.random() < .5 else -v) for v in rng.sample(range(1, nv + 1), min(3, nv))) + "\n")
    res = {"compile": (tc, mc, os.path.getsize(bd)), "nv": nv}
    for label, ba, da in QUERIES:
        if label not in labels:
            continue
        sub = lambda a: [x.replace("{c}", bd).replace("{q}", q) for x in a]
        b = timed([BDNNF] + sub(ba), tl)
        d = timed([DDNNIFE, "-i", nnf, "-t", str(nv)] + sub(da), tl)
        note = ""
        if label == "2-wise":
            for tool, r in (("bdnnf", b), ("ddnnife", d)):
                if r[0] is None:
                    continue
                f = f"{work}/{name}.{tool}.tw"
                open(f, "w").write(r[2])
                ok = subprocess.run([BDNNF, "twise-check", bd, f], capture_output=True).returncode == 0
                n = sum(1 for l in r[2].split("\n") if l.strip())
                note += f"{tool} {n}{'' if ok else ' (INCOMPLETE)'} "
        elif b[0] is not None and d[0] is not None and label in ("count", "features", "count queries",
                                                                  "SAT queries"):
            dl = [x for x in d[2].split("\n") if x.strip()]
            if label == "count":
                bo, do = b[2].split(), dl
            elif label == "features":
                bo = [x.split()[1] for x in b[2].split("\n") if x.strip() and not x.startswith("-")]
                do = [x.split(",")[1] for x in dl]
            else:
                bo = b[2].split()
                do = [x.rsplit(",", 1)[1] for x in dl]
                if label == "SAT queries":
                    do = ["SAT" if x == "true" else "UNSAT" for x in do]
            note = "agree" if bo == do else "DISAGREE"
        res[label] = (b[0], b[1], d[0], d[1], note.strip())
    return name, res


def fmt(t):
    return "TO" if t is None else (f"{t:.2f}s" if t < 10 else f"{t:.1f}s")


def main():
    ap = argparse.ArgumentParser(description="bdnnf vs ddnnife on shrunk competition CNFs")
    ap.add_argument("-j", "--jobs", type=int, default=11)
    ap.add_argument("--only", nargs="*", help="instance names")
    ap.add_argument("--timeout", type=int, default=120, help="per query, seconds")
    args = ap.parse_args()
    for b in (GANAK, BDNNF, DDNNIFE, CMS):
        if not os.access(b, os.X_OK):
            sys.exit(f"missing executable: {b}")
    insts = [i for i in INSTANCES if not args.only or i[0] in args.only]
    insts = [i for i in insts if os.path.exists(i[1]) or print(f"skipping {i[0]}: no {i[1]}")]
    work = tempfile.mkdtemp(prefix="bench_bdnnf_")
    t0 = time.time()
    with ProcessPoolExecutor(args.jobs) as pool:
        results = dict(pool.map(one, [(n, p, f, q, work, args.timeout) for n, p, f, q in insts]))
    print(f"\n{len(results)} instances, {time.time() - t0:.0f}s; work files in {work}\n")
    print("| instance | vars | compile | circuit |")
    print("|---|---|---|---|")
    for n, r in results.items():
        tc, mc, sz = r["compile"]
        print(f"| {n} | {r['nv']} | {fmt(tc)}, {mc or 0:.0f} MB | {sz / 1024:.0f} KB |")
    print()
    print("| query | instance | bdnnf | ddnnife | notes |")
    print("|---|---|---|---|---|")
    for label, _, _ in QUERIES:
        for n, r in results.items():
            if label not in r:
                continue
            bt, bm, dt, dm, note = r[label]
            print(f"| {label} | {n} | {fmt(bt)}, {bm or 0:.0f} MB | {fmt(dt)}, {dm or 0:.0f} MB | {note} |")
    print()
    print("| query | bdnnf total | ddnnife total | bdnnf faster | timeouts (bdnnf/ddnnife) |")
    print("|---|---|---|---|---|")
    for label, _, _ in QUERIES:
        rs = [r[label] for r in results.values() if label in r]
        both = [x for x in rs if x[0] is not None and x[2] is not None]
        print(f"| {label} | {fmt(sum(x[0] for x in both))} | {fmt(sum(x[2] for x in both))} | "
              f"{sum(x[0] < x[2] for x in both)}/{len(both)} | "
              f"{sum(x[0] is None for x in rs)}/{sum(x[2] is None for x in rs)} |")


if __name__ == "__main__":
    main()
