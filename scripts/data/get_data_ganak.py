#!/usr/bin/python3

import argparse
import collections
import csv
import decimal
import glob
import os
import re
import sqlite3 as sqlite_mod
import string
import sys

sys.set_int_max_str_digits(2000000)


##########################
# gpmc

#c o Components            = 339421
#c o conflicts             = 291050      (count 287680, sat 3370)
#c o decisions             = 643616      (count 335433, sat 308183)
def find_gpmc_time_cnt(fname):
    result = {}
    with open(fname, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith("c o Real time "):
                result["t"] = float(line.split()[5])
            elif line.startswith("c s exact arb int") or line.startswith("c s exact arb prec-sci"):
                val = line.split()[5]
                result["cnt"] = len(val) if len(val) > 1000 else decimal.Decimal(val).log10()
            elif line.startswith("c s log10-estimate"):
                result["mc_log10"] = float(line.split()[3])
            elif line.startswith("c o Components"):
                result["compsK"] = int(line.split()[4]) / 1000
            elif line.startswith("c o conflicts"):
                result["conflicts"] = int(line.split()[4])
            elif line.startswith("c o decisions"):
                result["decisionsK"] = int(line.split()[4]) / 1000
    return result


def approxmc_version(fname):
    aver = None
    cver = None
    with open(fname, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith("c ApproxMC SHA revision"):
                aver = line.split()[4]
            elif line.startswith("c CMS SHA revision"):
                cver = line.split()[4]
    if aver is not None:
        aver = "appmc" + aver[:8]
    if cver is not None:
        cver = cver[:8]
    return [aver, cver]


############################
## timeout wrapper
def timeout_parse(fname):
    t = None
    signal = None
    mem = None
    call = None
    solver = None
    page_faults = None
    with open(fname, "r") as f:
        for line in f:
            line = line.strip()
            if "Command terminated by signal" in line:
                signal = int(line.split()[4])
            elif "Minor (reclaiming a frame) page faults:" in line:
                page_faults = int(line.split()[6])
            elif "User time (seconds)" in line:
                t = float(line.split()[3])
            elif "Maximum resident set size (kbytes)" in line:
                mem = float(line.split()[5]) / 1000  # get it in MB
            elif "Command being timed" in line:
                call = " ".join(line.split()[3:])
                if "mc2022" in call:
                    call = call.split("mc2022_")[0]
                elif "mc2023" in call:
                    call = call.split("mc2023_")[0]
                else:
                    call = " ".join(call.split()[:-1])
                call = call.replace(" -t real", "")
                m = re.search(r"doalarm \d+", call)
                if m:
                    call = call.split(m.group(0))[1]
                else:
                    print("WARNING: no doalarm found in call: ", call)

                if "./ganak" in call:
                    solver = "ganak"
                if "./d4" in call:
                    solver = "d4"
                if "./approxmc" in call:
                    solver = "approxmc"
                if "./gpmc" in call:
                    solver = "gpmc"
                if "./gpmc-complex" in call:
                    solver = "gpmc"
                if "./KCBox" in call:
                    solver = "exactmc"
                if "./sharpSAT" in call:
                    solver = "sharptd"

                call = call.replace("././ganak ", "")
                call = call.replace("././d4-1d9cc6146f18b8 ", "")
                call = call.replace("././approxmc ", "")
                call = call.replace("././gpmc ", "")
                call = call.replace("././gpmc-complex ", "")
                call = call.replace("././KCBox-371eb601f2aa", "")
                call = call.replace("././sharpSAT", "")
                call = call.strip()

    if signal is None:
        signal = ""
    return {"t": t, "mem": mem, "call": call, "solver": solver,
            "page_faults": page_faults, "signal": signal}


# c o width 45
# c o CMD: timeout 60.000000s ./flow_cutter_pace17 ...
def sstd_treewidth(fname):
    result = {}
    with open(fname, "r") as f:
        for line in f:
            line = line.strip()
            if "c o width" in line:
                result["td_width"] = int(line.split()[3])
            elif "CMD: timeout" in line:
                result["td_time"] = float(line.split()[4][:-1])
    return result


def find_bad_solve(fname):
    mem_out = 0
    not_solved = True
    errored = 0
    with open(fname, "r") as f:
        for line in f:
            line = line.strip()
            if "bad_alloc" in line:
                mem_out = 1
            elif line.startswith("s SATISFIABLE"):
                not_solved = False
            elif line.startswith("s UNSATISFIABLE"):
                not_solved = False
            if "ERROR" in line or "assertion fail" in line.lower():
                errored = 1
    return mem_out, not_solved, errored


############################
## ganak — single-pass parser combining all per-line extractions
def parse_ganak_output(fname):
    result = {
        "errored": 0,
        "mem_out": 0,
        "not_solved": True,
        "cache_del_time": 0.0,
        "cubes_orig": 0,
        "cubes_final": 0,
    }
    aver = None
    cver = None

    with open(fname, "r") as f:
        for line in f:
            line = line.replace("[0m", "").replace("\x1b", "")
            line = ''.join(filter(lambda x: x in string.printable, line))
            line = line.strip()

            # Unconditional status checks
            if "ERROR" in line or "assertion fail" in line.lower():
                result["errored"] = 1
            if "Cannot allocate memory" in line or"bad_alloc" in line:
                result["mem_out"] = 1
            if line.startswith("s SATISFIABLE") or line.startswith("s UNSATISFIABLE"):
                result["not_solved"] = False

            # Mutually exclusive pattern matching
            if line.startswith("c o conflicts") and " :" in line:  # cryptominisat style
                result["conflicts"] = int(line.split()[4])
            elif line.startswith("c o conflicts"):
                result["conflicts"] = int(line.split()[3])
            elif line.startswith("c o decisions K"):
                result["decisionsK"] = int(line.split()[4])
            elif line.startswith("c GANAK SHA revision"):
                aver = line.split()[4]
            elif line.startswith("c CMS version"):
                cver = line.split()[3]
            elif line.startswith("c o GANAK SHA revision"):
                aver = line.split()[5]
            elif line.startswith("c p CMS version"):
                cver = line.split()[4]
            elif line.startswith("c o CMS revision"):
                cver = line.split()[4]
            elif line.startswith("c o buddy called /unsat ratio"):
                result["bdd_called"] = int(line.split()[6])
            elif line.startswith("c o buddy called"):
                result["bdd_called"] = int(line.split()[4])
            elif line.startswith("c o Sampling set size:") and "indepsz" not in result:
                indep_sz = line.split()[5].strip()
                indep_sz = 0 if indep_sz == "" else int(indep_sz)
                result["indepsz"] = 0 if indep_sz == 4294967295 else indep_sz
            elif line.startswith("c o opt ind size") and "newnvars" not in result:
                result["newnvars"] = int(line.split()[10])
            elif line.startswith("c o Opt sampling set size:") and "optindepsz" not in result:
                opt_indep_sz = line.split()[6].strip()
                opt_indep_sz = 0 if opt_indep_sz == "" else int(opt_indep_sz)
                result["optindepsz"] = 0 if opt_indep_sz == 4294967295 else opt_indep_sz
            elif line.startswith("c o CNF projection set size:"):
                result["origprojsz"] = int(line.split()[6])
            elif line.startswith("c o [extend-gates] Gates added to opt"):
                result["gates_extended"] = int(line.split()[8])
                result["gates_extend_t"] = float(line.split()[10])
            elif line.startswith("c o [arjun-extend] Extend finished"):
                orig = int(line.split()[7])
                final = int(line.split()[10])
                result["padoa_extended"] = final - orig
                result["padoa_extend_t"] = float(line.split()[14])
            elif line.startswith("c o [arjun] Start unknown size:"):
                result["unknsz"] = int(line.split()[6])
            elif "c o [arjun] backward round finished" in line:
                result["backwtime"] = float(line.split()[10])
            elif line.startswith("c o  ") and "% total" in line:
                result["backboneT"] = result.get("backboneT", 0) + float(line.split()[2])
            elif line.startswith("c o Arjun T:"):
                result["arjuntime"] = float(line.split()[4])
            elif line.startswith("c o [td] iter") and "best bag" in line and "td_width" not in result:
                result["td_width"] = int(line.split()[7]) - 1
                result["td_time"] = float(line.split()[12])
            elif line.startswith("c o [td] iter") and "width:" in line and "td_width" not in result:
                result["td_width"] = int(line.split()[6]) - 1
                result["td_time"] = float(line.split()[11])
            elif line.startswith("c o [td] Primal graph"):
                result["primal_density"] = float(line.split()[10])
                result["primal_edge_var_ratio"] = float(line.split()[12])
            elif line.startswith("c o cache miss rate"):
                result["cache_miss_rate"] = float(line.split()[5])
            elif line.startswith("c o cache K (lookup/ stores/ hits/ dels)"):
                result["compsK"] = float(line.split()[8])
            elif line.startswith("c o avg hit/store num vars"):
                result["cache_avg_hit_vars"] = float(line.split()[6])
                result["cache_avg_store_vars"] = float(line.split()[8])
            elif line.startswith("c o deletion done. T:"):
                result["cache_del_time"] += float(line.split()[5])
            elif line.startswith("c o cubes orig:"):
                result["cubes_orig"] += int(line.split()[4])
                result["cubes_final"] += int(line.split()[6])
            elif line.startswith("c o Num restarts:"):
                result["restarts"] = int(line.split()[4])
            elif line.startswith("c o [rst-cube] Num restarts:"):
                result["restarts"] = int(line.split()[5])
            elif line.startswith("c o sat called/sat/unsat/conflK"):
                result["sat_called"] = int(line.split()[4])
            elif line.startswith("c o sat call/sat/unsat/confl/rst"):
                result["sat_called"] = int(line.split()[4])
                result["satrst"] = int(line.split()[8])
            elif line.startswith("c o Bin irred/red") and "irred_bin" not in result:
                result["irred_bin"] = int(line.split()[4])
            elif line.startswith("c o Long irred cls/tri") and "irred_long" not in result:
                result["irred_long"] = int(line.split()[5])
                result["irred_tri"] = int(line.split()[6])
            elif line.startswith("c s log10-estimate"):
                result["mc_log10"] = float(line.split()[3])

    if aver is not None:
        aver = aver[:8]
    if cver is not None:
        cver = cver[:8]
    result["solverver"] = ["ganak", "%s-%s" % (aver, cver)]

    irred_bin  = result.get("irred_bin",  0)
    irred_long = result.get("irred_long", 0)
    irred_tri  = result.get("irred_tri",  0)
    if irred_bin or irred_long or irred_tri:
        result["irred_cls"] = irred_bin + irred_long + irred_tri

    return result


_SIMP_STATS_RE = re.compile(
    r'c o \[simp-stats\] (BEFORE|AFTER) (\S+) '
    r'irred_bins (\d+) irred_long_cls (\d+) irred_long_lits (\d+) free_vars (\d+)'
    r'(?:\s+elimed_vars (\d+) replaced_vars (\d+) units (\d+) mem_MB (\d+) T:\s*([\d.]+))?'
    r'(?:\s+depth (\d+))?'
)


def parse_simp_stats(fname):
    """Parse [simp-stats] BEFORE/AFTER pairs from a ganak output file.
    Returns a list of dicts, one per matched BEFORE/AFTER pair.

    Each record carries `depth` (nesting level of the wrapper at emission
    time). Top-level steps are depth=0; nested children are depth>=1 and
    their work is already counted by the parent — sum only depth=0 for
    total-time accounting.

    Older logs without the `depth N` suffix get depth=0 (treated as flat).
    """
    # pending[(name, depth)] -> stack of pending BEFOREs awaiting matching AFTER.
    # LIFO matching is required when the same-name step nests (e.g. clean-cls
    # called from inside another clean-cls dispatch).
    pending = {}
    step_counts = {}  # name -> number of completed pairs so far
    records = []
    last_completed_step = "none"  # tracks prev_step

    with open(fname, "r") as f:
        for line in f:
            line = line.replace("[0m", "").replace("\x1b", "").strip()
            m = _SIMP_STATS_RE.search(line)
            if not m:
                continue
            kind = m.group(1)
            name = m.group(2)
            # Core fields (always present)
            core = (int(m.group(3)), int(m.group(4)), int(m.group(5)), int(m.group(6)))
            # Extended fields (present in newer logs)
            ext = None
            if m.group(7) is not None:
                ext = (int(m.group(7)), int(m.group(8)), int(m.group(9)),
                       int(m.group(10)), float(m.group(11)))
            depth = int(m.group(12)) if m.group(12) is not None else 0

            key = (name, depth)
            if kind == "BEFORE":
                pending.setdefault(key, []).append((core, ext))
            elif kind == "AFTER":
                if pending.get(key):
                    before_core, before_ext = pending[key].pop()
                    occurrence = step_counts.get(name, 0)
                    step_counts[name] = occurrence + 1
                    step_num = len(records)
                    rec = {
                        "step_num":  step_num,
                        "occurrence": occurrence,
                        "name":      name,
                        "depth":     depth,
                        "prev_step": last_completed_step,
                        "before_irred_bins":      before_core[0],
                        "before_irred_long_cls":  before_core[1],
                        "before_irred_long_lits": before_core[2],
                        "before_free_vars":       before_core[3],
                        "after_irred_bins":       core[0],
                        "after_irred_long_cls":   core[1],
                        "after_irred_long_lits":  core[2],
                        "after_free_vars":        core[3],
                        "delta_irred_bins":       core[0] - before_core[0],
                        "delta_irred_long_cls":   core[1] - before_core[1],
                        "delta_irred_long_lits":  core[2] - before_core[2],
                        "delta_free_vars":        core[3] - before_core[3],
                    }
                    # Extended fields
                    if ext is not None and before_ext is not None:
                        rec["before_elimed_vars"]   = before_ext[0]
                        rec["before_replaced_vars"] = before_ext[1]
                        rec["before_units"]         = before_ext[2]
                        rec["before_mem_MB"]        = before_ext[3]
                        rec["after_elimed_vars"]    = ext[0]
                        rec["after_replaced_vars"]  = ext[1]
                        rec["after_units"]          = ext[2]
                        rec["after_mem_MB"]         = ext[3]
                        rec["delta_elimed_vars"]    = ext[0] - before_ext[0]
                        rec["delta_replaced_vars"]  = ext[1] - before_ext[1]
                        rec["delta_units"]          = ext[2] - before_ext[2]
                        rec["delta_mem_MB"]         = ext[3] - before_ext[3]
                        rec["step_time"]            = ext[4] - before_ext[4]
                    records.append(rec)
                    last_completed_step = name
    return records


def load_already_parsed(db_path):
    """Return set of (dirname, fname) tuples already in data.sqlite3."""
    if not os.path.exists(db_path):
        return set()
    try:
        conn = sqlite_mod.connect(db_path)
        cur = conn.cursor()
        cur.execute("SELECT dirname, fname FROM data")
        result = set(cur.fetchall())
        conn.close()
        return result
    except Exception:
        return set()


def main():
    parser = argparse.ArgumentParser(description="Parse ganak output files into CSV")
    parser.add_argument("--files", default="out-ganak*/*cnf*",
                        help="Glob pattern for input files (default: 'out-ganak*/*cnf*')")
    parser.add_argument("--verbose", action="store_true",
                        help="Print each file being parsed")
    parser.add_argument("--reparse", action="store_true",
                        help="Delete data.sqlite3 and reparse all files from scratch")
    args = parser.parse_args()

    if args.reparse and os.path.exists("data.sqlite3"):
        os.remove("data.sqlite3")
        print("Deleted data.sqlite3, reparsing all files")

    already_parsed = load_already_parsed("data.sqlite3")
    if already_parsed:
        skip_counts = collections.Counter(dirname for dirname, _ in already_parsed)
        print(f"Skipping {len(already_parsed)} already-parsed files from data.sqlite3:")
        for d in sorted(skip_counts):
            print(f"  {d}: {skip_counts[d]} files")

    file_list = glob.glob(args.files)
    files = {}
    for f in file_list:
        if ".csv" in f:
            continue

        dirname = f.split("/")[0]
        if "competitors" in dirname:
            continue
        full_fname = f.split("/")[1].split(".cnf")[0] + ".cnf"
        base = dirname + "/" + full_fname

        if (dirname, full_fname) in already_parsed:
            if args.verbose:
                print(f"skipping (already in DB): {base}")
            continue

        if args.verbose:
            print("parsing file: ", f)

        if base not in files:
            files[base] = {}
        files[base]["dirname"] = dirname
        files[base]["fname"] = full_fname

        if f.endswith(".timeout") or ".timeout_" in f:
            tp = timeout_parse(f)
            files[base]["timeout_t"] = tp["t"]
            files[base]["timeout_mem"] = tp["mem"]
            files[base]["timeout_call"] = tp["call"]
            files[base]["page_faults"] = tp["page_faults"]
            files[base]["signal"] = tp["signal"]
            if "solver" not in files[base] or files[base]["solver"] is None:
                files[base]["solver"] = tp["solver"]
            continue

        if f.endswith(".out_ganak") or f.endswith(".out_arjun") or f.endswith(".out"):
            files[base]["solver"] = "ganak"
            files[base].update(parse_ganak_output(f))
            simp = parse_simp_stats(f)
            if simp:
                files[base]["simp_stats"] = simp
        if ".out_approxmc" in f:
            files[base]["solver"] = "approxmc"
            files[base]["solverver"] = approxmc_version(f)
        if ".out_gpmc" in f:
            files[base]["solver"] = "gpmc"
            gpmc = find_gpmc_time_cnt(f)
            files[base]["conflicts"] = gpmc.get("conflicts")
            files[base]["compsK"] = gpmc.get("compsK")
            files[base]["decisionsK"] = gpmc.get("decisionsK")
            files[base]["solverver"] = ["gpmc", "gpmc"]
        if ".out_sharptd" in f:
            files[base]["solver"] = "sharptd"
            files[base]["solverver"] = ["sharptd", "sharptd"]
            td = sstd_treewidth(f)
            files[base]["td_width"] = td.get("td_width")
            files[base]["td_time"] = td.get("td_time")
        if ".out_d4" in f:
            files[base]["solver"] = "d4"
            files[base]["solverver"] = ["d4", "d4"]
        if ".out_dsharp" in f:
            files[base]["solver"] = "dsharp"
            files[base]["solverver"] = ["dsharp", "dsharp"]
        if ".out_minic2d" in f:
            files[base]["solver"] = "minic2d"
            files[base]["solverver"] = ["minic2d", "minic2d"]
        if ".out_exactmc" in f:
            files[base]["solver"] = "exactmc"
            files[base]["solverver"] = ["exactmc", "exactmc"]

        if ".out" in f:
            mem_out, not_solved, errored = find_bad_solve(f)
            files[base]["mem_out"] = max(files[base].get("mem_out", 0), mem_out)
            files[base]["not_solved"] = not_solved
            files[base]["errored"] = max(files[base].get("errored", 0), errored)

    if not files:
        print("No new files to insert into data.sqlite3")
        return

    new_counts = collections.Counter(v["dirname"] for v in files.values())
    print(f"Parsing {len(files)} new files:")
    for d in sorted(new_counts):
        print(f"  {d}: {new_counts[d]} files")

    cols = ["solver", "dirname", "fname", "mem_out", "errored", "ganak_time", "ganak_mem_MB",
            "ganak_call", "page_faults", "signal", "ganak_ver", "conflicts", "decisionsK",
            "compsK", "primal_density", "primal_edge_var_ratio", "td_width", "td_time",
            "arjun_time", "backboneT", "backwardT", "indepsz", "optindepsz", "origprojsz",
            "new_nvars", "unknsz", "cache_del_time", "cache_avg_hit_vars",
            "cache_avg_store_vars", "cache_miss_rate", "bdd_called", "sat_called",
            "sat_rst", "rst", "cubes_orig", "cubes_final", "gates_extended",
            "gates_extend_t", "padoa_extended", "padoa_extend_t", "timeout_t",
            "irred_bin", "irred_long", "irred_tri", "irred_cls"]

    def g(d, key):
        v = d.get(key)
        return "" if v is None else v

    with open("mydata.csv", "w", newline="") as out:
        writer = csv.writer(out)
        writer.writerow(cols)
        for _, f in files.items():
            if "not_solved" not in f:
                print("WARNING no 'out' file parsed for, skipping: ", f["fname"])
                continue
            if "solver" not in f:
                print("oops, solver not found, that's wrong")
                print(f)
                exit(-1)
            if "timeout_t" not in f:
                print("timeout not parsed for f: ", f)
                exit(-1)

            ganak_time = "" if (f["timeout_t"] is None or f["not_solved"]) else f["timeout_t"]

            solverver = ""
            if "solverver" in f and f["solverver"] != [None, None]:
                solverver = "%s-%s" % (f["solverver"][0], f["solverver"][1])

            writer.writerow([
                g(f, "solver"),
                g(f, "dirname"),
                g(f, "fname"),
                g(f, "mem_out"),
                g(f, "errored"),
                ganak_time,
                g(f, "timeout_mem"),
                g(f, "timeout_call"),
                g(f, "page_faults"),
                g(f, "signal"),
                solverver,
                g(f, "conflicts"),
                g(f, "decisionsK"),
                g(f, "compsK"),
                g(f, "primal_density"),
                g(f, "primal_edge_var_ratio"),
                g(f, "td_width"),
                g(f, "td_time"),
                g(f, "arjuntime"),
                g(f, "backboneT"),
                g(f, "backwtime"),
                g(f, "indepsz"),
                g(f, "optindepsz"),
                g(f, "origprojsz"),
                g(f, "newnvars"),
                g(f, "unknsz"),
                g(f, "cache_del_time"),
                g(f, "cache_avg_hit_vars"),
                g(f, "cache_avg_store_vars"),
                g(f, "cache_miss_rate"),
                g(f, "bdd_called"),
                g(f, "sat_called"),
                g(f, "satrst"),
                g(f, "restarts"),
                g(f, "cubes_orig"),
                g(f, "cubes_final"),
                g(f, "gates_extended"),
                g(f, "gates_extend_t"),
                g(f, "padoa_extended"),
                g(f, "padoa_extend_t"),
                g(f, "timeout_t"),
                g(f, "irred_bin"),
                g(f, "irred_long"),
                g(f, "irred_tri"),
                g(f, "irred_cls"),
            ])

    preproc_cols = [
        "dirname", "fname", "step_num", "occurrence", "name", "depth", "prev_step",
        "before_irred_bins", "before_irred_long_cls", "before_irred_long_lits", "before_free_vars",
        "after_irred_bins", "after_irred_long_cls", "after_irred_long_lits", "after_free_vars",
        "delta_irred_bins", "delta_irred_long_cls", "delta_irred_long_lits", "delta_free_vars",
        # Extended fields (newer log format)
        "before_elimed_vars", "before_replaced_vars", "before_units", "before_mem_MB",
        "after_elimed_vars", "after_replaced_vars", "after_units", "after_mem_MB",
        "delta_elimed_vars", "delta_replaced_vars", "delta_units", "delta_mem_MB",
        "step_time",
    ]

    def g_rec(rec, key):
        v = rec.get(key)
        return "" if v is None else v

    with open("preproc.csv", "w", newline="") as out:
        writer = csv.writer(out)
        writer.writerow(preproc_cols)
        for _, f in files.items():
            if "simp_stats" not in f:
                continue
            for rec in f["simp_stats"]:
                writer.writerow([
                    f["dirname"],
                    f["fname"],
                    rec["step_num"],
                    rec["occurrence"],
                    rec["name"],
                    rec.get("depth", 0),
                    rec.get("prev_step", "none"),
                    rec["before_irred_bins"],
                    rec["before_irred_long_cls"],
                    rec["before_irred_long_lits"],
                    rec["before_free_vars"],
                    rec["after_irred_bins"],
                    rec["after_irred_long_cls"],
                    rec["after_irred_long_lits"],
                    rec["after_free_vars"],
                    rec["delta_irred_bins"],
                    rec["delta_irred_long_cls"],
                    rec["delta_irred_long_lits"],
                    rec["delta_free_vars"],
                    g_rec(rec, "before_elimed_vars"),
                    g_rec(rec, "before_replaced_vars"),
                    g_rec(rec, "before_units"),
                    g_rec(rec, "before_mem_MB"),
                    g_rec(rec, "after_elimed_vars"),
                    g_rec(rec, "after_replaced_vars"),
                    g_rec(rec, "after_units"),
                    g_rec(rec, "after_mem_MB"),
                    g_rec(rec, "delta_elimed_vars"),
                    g_rec(rec, "delta_replaced_vars"),
                    g_rec(rec, "delta_units"),
                    g_rec(rec, "delta_mem_MB"),
                    g_rec(rec, "step_time"),
                ])

    def n(v):
        """Convert empty string to None for sqlite NULL."""
        return None if v == "" else v

    conn = sqlite_mod.connect("data.sqlite3")
    conn.executescript("""
        CREATE TABLE IF NOT EXISTS data (
          solver STRING NOT NULL,
          dirname STRING NOT NULL,
          fname STRING NOT NULL,
          mem_out INT,
          errored INT,
          ganak_time FLOAT,
          ganak_mem_MB FLOAT,
          ganak_call TEXT NOT NULL,
          page_faults INT NOT NULL,
          signal INT,
          ganak_ver TEXT NOT NULL,
          conflicts INT,
          decisionsK INT,
          compsK INT,
          primal_density FLOAT,
          primal_edge_var_ratio FLOAT,
          td_width INT,
          td_time FLOAT,
          arjun_time FLOAT,
          backbone_time FLOAT,
          backward_time FLOAT,
          indep_sz INT,
          opt_indep_sz INT,
          orig_proj_sz INT,
          new_nvars INT,
          unkn_sz INT,
          cache_del_time FLOAT,
          cache_avg_hit_vars FLOAT,
          cache_avg_store_vars FLOAT,
          cache_miss_rate FLOAT,
          bdd_called INT,
          sat_called INT,
          sat_rst INT,
          restarts INT,
          cubes_orig INT,
          cubes_final INT,
          gates_extended INT,
          gates_extend_t INT,
          padoa_extended INT,
          padoa_extend_t INT,
          timeout_t FLOAT,
          irred_bin INT,
          irred_long INT,
          irred_tri INT,
          irred_cls INT,
          mc_log10 FLOAT
        );
        CREATE TABLE IF NOT EXISTS preproc (
          dirname STRING NOT NULL,
          fname STRING NOT NULL,
          step_num INT NOT NULL,
          occurrence INT NOT NULL,
          name STRING NOT NULL,
          depth INT NOT NULL DEFAULT 0,
          prev_step STRING NOT NULL DEFAULT 'none',
          before_irred_bins INT,
          before_irred_long_cls INT,
          before_irred_long_lits INT,
          before_free_vars INT,
          after_irred_bins INT,
          after_irred_long_cls INT,
          after_irred_long_lits INT,
          after_free_vars INT,
          delta_irred_bins INT,
          delta_irred_long_cls INT,
          delta_irred_long_lits INT,
          delta_free_vars INT,
          before_elimed_vars INT,
          before_replaced_vars INT,
          before_units INT,
          before_mem_MB INT,
          after_elimed_vars INT,
          after_replaced_vars INT,
          after_units INT,
          after_mem_MB INT,
          delta_elimed_vars INT,
          delta_replaced_vars INT,
          delta_units INT,
          delta_mem_MB INT,
          step_time FLOAT
        );
    """)

    data_rows = []
    for _, f in files.items():
        if "not_solved" not in f:
            print("WARNING no 'out' file parsed for, skipping: ", f["fname"])
            continue
        if "solver" not in f:
            print("oops, solver not found, that's wrong")
            print(f)
            exit(-1)
        if "timeout_t" not in f:
            print("timeout not parsed for f: ", f)
            exit(-1)

        ganak_time = None if (f["timeout_t"] is None or f["not_solved"]) else f["timeout_t"]

        solverver = ""
        if "solverver" in f and f["solverver"] != [None, None]:
            solverver = "%s-%s" % (f["solverver"][0], f["solverver"][1])

        data_rows.append((
            f.get("solver", ""),
            f.get("dirname", ""),
            f.get("fname", ""),
            n(f.get("mem_out", "")),
            n(f.get("errored", "")),
            ganak_time,
            n(f.get("timeout_mem", "")),
            f.get("timeout_call", "") or "",
            f.get("page_faults", "") or "",
            n(f.get("signal", "")),
            solverver,
            n(f.get("conflicts", "")),
            n(f.get("decisionsK", "")),
            n(f.get("compsK", "")),
            n(f.get("primal_density", "")),
            n(f.get("primal_edge_var_ratio", "")),
            n(f.get("td_width", "")),
            n(f.get("td_time", "")),
            n(f.get("arjuntime", "")),
            n(f.get("backboneT", "")),
            n(f.get("backwtime", "")),
            n(f.get("indepsz", "")),
            n(f.get("optindepsz", "")),
            n(f.get("origprojsz", "")),
            n(f.get("newnvars", "")),
            n(f.get("unknsz", "")),
            n(f.get("cache_del_time", "")),
            n(f.get("cache_avg_hit_vars", "")),
            n(f.get("cache_avg_store_vars", "")),
            n(f.get("cache_miss_rate", "")),
            n(f.get("bdd_called", "")),
            n(f.get("sat_called", "")),
            n(f.get("satrst", "")),
            n(f.get("restarts", "")),
            n(f.get("cubes_orig", "")),
            n(f.get("cubes_final", "")),
            n(f.get("gates_extended", "")),
            n(f.get("gates_extend_t", "")),
            n(f.get("padoa_extended", "")),
            n(f.get("padoa_extend_t", "")),
            n(f.get("timeout_t", "")),
            n(f.get("irred_bin", "")),
            n(f.get("irred_long", "")),
            n(f.get("irred_tri", "")),
            n(f.get("irred_cls", "")),
            n(f.get("mc_log10", "")),
        ))

    conn.executemany(
        "INSERT INTO data VALUES (" + ",".join(["?"] * 46) + ")",
        data_rows
    )

    preproc_rows = []
    for _, f in files.items():
        if "simp_stats" not in f:
            continue
        for rec in f["simp_stats"]:
            preproc_rows.append((
                f["dirname"],
                f["fname"],
                rec["step_num"],
                rec["occurrence"],
                rec["name"],
                rec.get("depth", 0),
                rec.get("prev_step", "none"),
                n(rec.get("before_irred_bins")),
                n(rec.get("before_irred_long_cls")),
                n(rec.get("before_irred_long_lits")),
                n(rec.get("before_free_vars")),
                n(rec.get("after_irred_bins")),
                n(rec.get("after_irred_long_cls")),
                n(rec.get("after_irred_long_lits")),
                n(rec.get("after_free_vars")),
                n(rec.get("delta_irred_bins")),
                n(rec.get("delta_irred_long_cls")),
                n(rec.get("delta_irred_long_lits")),
                n(rec.get("delta_free_vars")),
                n(rec.get("before_elimed_vars")),
                n(rec.get("before_replaced_vars")),
                n(rec.get("before_units")),
                n(rec.get("before_mem_MB")),
                n(rec.get("after_elimed_vars")),
                n(rec.get("after_replaced_vars")),
                n(rec.get("after_units")),
                n(rec.get("after_mem_MB")),
                n(rec.get("delta_elimed_vars")),
                n(rec.get("delta_replaced_vars")),
                n(rec.get("delta_units")),
                n(rec.get("delta_mem_MB")),
                n(rec.get("step_time")),
            ))

    conn.executemany(
        "INSERT INTO preproc VALUES (" + ",".join(["?"] * 32) + ")",
        preproc_rows
    )

    conn.commit()
    conn.close()
    print(f"Inserted {len(data_rows)} new rows into data, {len(preproc_rows)} rows into preproc")


if __name__ == "__main__":
    main()
