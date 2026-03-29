#!/usr/bin/python3

import argparse
import csv
import decimal
import glob
import os
import string
import subprocess
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
                if "doalarm 3600" in call:
                    call = call.split("doalarm 3600")[1]
                elif "doalarm 60" in call:
                    call = call.split("doalarm 60")[1]
                else:
                    call = call.split("doalarm 900")[1]

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
            if "bad_alloc" in line:
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
            elif line.startswith("c o Sampling set size:"):
                indep_sz = line.split()[5].strip()
                indep_sz = 0 if indep_sz == "" else int(indep_sz)
                result["indepsz"] = 0 if indep_sz == 4294967295 else indep_sz
            elif line.startswith("c o opt ind size"):
                result["newnvars"] = int(line.split()[10])
            elif line.startswith("c o Opt sampling set size:"):
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
            elif line.startswith("c o [td] iter") and "best bag" in line:
                result["td_width"] = int(line.split()[7]) - 1
                result["td_time"] = float(line.split()[12])
            elif line.startswith("c o [td] iter") and "width:" in line:
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

    if aver is not None:
        aver = aver[:8]
    if cver is not None:
        cver = cver[:8]
    result["solverver"] = ["ganak", "%s-%s" % (aver, cver)]

    return result


def main():
    parser = argparse.ArgumentParser(description="Parse ganak output files into CSV")
    parser.add_argument("--files", default="out-ganak*/*cnf*",
                        help="Glob pattern for input files (default: 'out-ganak*/*cnf*')")
    parser.add_argument("--verbose", action="store_true",
                        help="Print each file being parsed")
    args = parser.parse_args()

    file_list = glob.glob(args.files)
    files = {}
    for f in file_list:
        if ".csv" in f:
            continue
        if args.verbose:
            print("parsing file: ", f)

        dirname = f.split("/")[0]
        if "competitors" in dirname:
            continue
        full_fname = f.split("/")[1].split(".cnf")[0] + ".cnf"
        base = dirname + "/" + full_fname
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

        if f.endswith(".out_ganak") or f.endswith(".out"):
            files[base]["solver"] = "ganak"
            files[base].update(parse_ganak_output(f))
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

    cols = ["solver", "dirname", "fname", "mem_out", "errored", "ganak_time", "ganak_mem_MB",
            "ganak_call", "page_faults", "signal", "ganak_ver", "conflicts", "decisionsK",
            "compsK", "primal_density", "primal_edge_var_ratio", "td_width", "td_time",
            "arjun_time", "backboneT", "backwardT", "indepsz", "optindepsz", "origprojsz",
            "new_nvars", "unknsz", "cache_del_time", "cache_avg_hit_vars",
            "cache_avg_store_vars", "cache_miss_rate", "bdd_called", "sat_called",
            "sat_rst", "rst", "cubes_orig", "cubes_final", "gates_extended",
            "gates_extend_t", "padoa_extended", "padoa_extend_t", "timeout_t"]

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
            ])

    with open("ganak.sqlite") as sql_f:
        subprocess.run(["sqlite3", "data.sqlite3"], stdin=sql_f, check=True)


if __name__ == "__main__":
    main()
