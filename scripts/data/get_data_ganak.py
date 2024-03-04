#!/usr/bin/python3

import os
import glob
import decimal
import sys

sys.set_int_max_str_digits(2000000)

def lit_cnt(fname):
    num_cls = None
    num_bin_cls = None
    lits = None
    vars = None
    with open(fname, "r") as f:
        for line in f:
            line=line.strip()
            # print(line)
            # print(line.split())
            if "num cls" in line:
                num_cls = int(line.split()[2])
            if "num bin cls" in line:
                num_bin_cls = int(line.split()[3])
            if "num (non-unit) lits" in line:
                lits = int(line.split()[3])
            if "num (non-set) vars" in line:
                vars = int(line.split()[3])

    return [vars, num_cls, num_bin_cls, lits]

##########################
# appmx, etc.

def find_approxmc_time_cnt(fname):
    t = None
    cnt = None
    with open(fname, "r") as f:
        for line in f:
            line = line.strip()
            if "c [appmc+arjun] Total time" in line:
                t = float(line.split()[4])
            if "s mc" in line:
                if len(line.split()[2]) > 1000 : cnt = len(line.split()[2])
                else: cnt = decimal.Decimal(line.split()[2]).log10()

    return [t,cnt]

def find_exactmc_time_cnt(fname):
    t = None
    cnt = None
    with open(fname, "r") as f:
        for line in f:
            line = line.strip()
            if "Total time cost" in line:
                t = float(line.split()[3])
            if "Number of models" in line:
                if len(line.split()[3]) > 1000: cnt = len(line.split()[3])
                cnt = decimal.Decimal(line.split()[3]).log10()

    return [t,cnt]

def find_minic2d_time_cnt(fname):
    t = None
    cnt = None
    with open(fname, "r") as f:
        for line in f:
            line = line.strip()
            if "Total Time" in line:
                t = float(line.split()[2].strip("s"))
            if "Counting..." in line:
                if len(line.split()[1]) > 1000: cnt = len(line.split()[1])
                else: cnt = decimal.Decimal(line.split()[1]).log10()

    return [t,cnt]

def find_d4_time_cnt(fname):
    t = None
    cnt = None
    with open(fname, "r") as f:
        for line in f:
            line = line.strip()
            if "c Final time:" in line:
                t = float(line.split()[3])
            if len(line) >=3 and line[:2] == "s ":
                if len(line.split()[1]) > 1000: cnt = line.split()[1]
                else: cnt = decimal.Decimal(line.split()[1]).log10()

    # print("t:", t, "cnt: ", cnt)
    return [t,cnt]

def find_gpmc_time_cnt(fname):
    t = None
    cnt = None
    with open(fname, "r") as f:
        for line in f:
            line = line.strip()
            if "c o Real time " in line:
                t = float(line.split()[5])
            if "c s exact arb int" in line:
                if len(line.split()[5]) > 1000: cnt = len(line.split()[5])
                else: cnt = decimal.Decimal(line.split()[5]).log10()

    return [t,cnt]

def find_dsharp_time_cnt(fname):
    t = None
    cnt = None
    with open(fname, "r") as f:
        for line in f:
            line = line.strip()
            if "Runtime: " in line:
                t = float(line.split(":")[1].strip("s"))
            if "#SAT (full)" in line:
                if len(line.split()[2]) > 1000: cnt = len(line.split()[2])
                else: cnt = decimal.Decimal(line.split()[2]).log10()

    return [t,cnt]

# when only preprocessing finds an answer, we don't have "c o Solved"
# Instead, we have only "c o Preprocessed. 469.010186778s Vars: 0 Clauses: 0 Free vars: 1"
def find_sharpsat_time_cnt(fname):
    t = None
    cnt = None
    prepro_t = None
    with open(fname, "r") as f:
        for line in f:
            line = line.strip()
            if "c o Solved." in line:
                t = float(line.split()[3])
            if "c o Preprocessed." in line:
                d = line.split()[3]
                d = d.strip('s')
                prepro_t = float(d)
            if "c s exact arb int" in line:
                if len(line.split()[5]) > 1000: cnt = len(line.split()[5])
                else: cnt = decimal.Decimal(line.split()[5]).log10()

    if cnt is not None and t is None:
      t = prepro_t
    return [t,cnt]

def approxmc_version(fname):
    aver = None
    cver = None
    with open(fname, "r") as f:
        for line in f:
            line = line.strip()
            if "c ApproxMC SHA revision" in line:
                aver = line.split()[4]
            if "c CMS SHA revision" in line:
                cver = line.split()[4]

    if aver is not None:
        aver = "appmc"+aver[:8]
    if cver is not None:
        cver = cver[:8]
    return [aver, cver]

############################
## ganak

def find_ganak_time_cnt(fname):
    t = None
    cnt = None
    with open(fname, "r") as f:
        for line in f:
            line = line.strip()
            if "c Time" in line:
                t = float(line.split()[2])
            if "c o Total time [Arjun+GANAK]:" in line:
                t = float(line.split()[5])
            if "s mc" in line:
                cnt = decimal.Decimal(line.split()[2])

    return [t,cnt]


def timeout_parse(fname):
    t = None
    m = None
    call = None
    with open(fname, "r") as f:
        for line in f:
            line = line.strip()
            if "User time (seconds)" in line:
                t = float(line.split()[3])
            if "Maximum resident set size (kbytes)" in line:
                m = float(line.split()[5])/(1000) # get it in MB
            if "Command being timed" in line:
                call= " ".join(line.split()[3:])
                if "mc2022" in call:
                    call = call.split("mc2022_")[0]
                elif "mc2023" in call:
                    call = call.split("mc2023_")[0]
                else:
                    assert False
                call = call.replace(" -t real", "")
                call = call.split("doalarm 3600")[1]
                call = call.replace("././ganak ", "")
                call = call.replace("././d4-1d9cc6146f18b8 ", "")
                call = call.replace("././approxmc ", "")
                call = call.replace("././gpmc ", "")
                call = call.replace("././KCBox-371eb601f2aa", "")
                call = call.replace("././sharpSAT", "")
                call = call.strip()

    return [t, m, call]

# c o conflicts                      2503077      -- confl/s:    1434.45
def ganak_conflicts(fname) -> int:
    confl = 0
    with open(fname, "r") as f:
        for line in f:
            line = line.strip()
            if "c o conflicts   " in line:
                confl = int(line.split()[3])

    return confl


# c o cache K (lookup/ stores/ hits) 67229  34594  32634   -- Klookup/s:  38.53
def ganak_comps(fname) -> int:
    comps = 0
    with open(fname, "r") as f:
        for line in f:
            line = line.strip()
            if "c o cache K (lookup/ stores/ hits)" in line:
                comps = int(line.split()[7])

    return comps


def ganak_version(fname):
    aver = None
    cver = None
    with open(fname, "r") as f:
        for line in f:
            line = line.strip()
            if "c GANAK SHA revision" in line:
                aver = line.split()[4]
            if "c CMS version" in line:
                cver = line.split()[3]
            if "c o GANAK SHA revision" in line:
                aver = line.split()[5]
            if "c p CMS version" in line:
                cver = line.split()[4]

    if aver == "74816f58aa522ec39ed7a9dc118b9abadfb68f09":
        aver = "f6789ccc62c9748f03198467d2d24d2136901b1b"

    if aver is not None:
        aver = aver[:8]
    if cver is not None:
        cver = cver[:8]
    return ["ganak", "%s-%s" % (aver,cver)]


file_list = glob.glob("out-ganak-*/*cnf*")
files = {}
for f in file_list:
    if ".csv" in f:
        continue
    dirname = f.split("/")[0]
    if "competitors" in dirname:
        continue
    full_fname = f.split("/")[1].split(".cnf")[0] + ".cnf"
    base = dirname+"/"+full_fname
    if base not in files:
        files[base] = {}
    files[base]["dirname"] = dirname
    files[base]["fname"] = full_fname
    # print("Dealing with dir: %s fname: %s" % (dirname, full_fname))

    if ".timeout_" in f:
        files[base]["solvertout"] = timeout_parse(f)

    if ".out_ganak" in f:
        files[base]["solver"] = "ganak"
        files[base]["solvertime"] = find_ganak_time_cnt(f)
        files[base]["solverver"] = ganak_version(f)
        files[base]["confls"] = ganak_conflicts(f)
        files[base]["comps"] = ganak_comps(f)
    if ".out_approxmc" in f:
        files[base]["solver"] = "approxmc"
        files[base]["solvertime"] = find_approxmc_time_cnt(f)
        files[base]["solverver"] = approxmc_version(f)
        files[base]["confls"] = 0
        files[base]["comps"] = 0
    if ".out_gpmc" in f:
        files[base]["solver"] = "gpmc"
        files[base]["solvertime"] = find_gpmc_time_cnt(f)
        files[base]["solverver"] = ["gpmc", "gpmc"]
        files[base]["confls"] = 0
        files[base]["comps"] = 0
    if ".out_sharptd" in f:
        files[base]["solver"] = "sharptd"
        files[base]["solvertime"] = find_sharpsat_time_cnt(f)
        files[base]["solverver"] = ["sharptd", "sharptd"]
        files[base]["confls"] = 0
        files[base]["comps"] = 0
    if ".out_d4" in f:
        files[base]["solver"] = "d4"
        files[base]["solvertime"] = find_d4_time_cnt(f)
        files[base]["solverver"] = ["d4", "d4"]
        files[base]["confls"] = 0
        files[base]["comps"] = 0
    if ".out_dsharp" in f:
        files[base]["solver"] = "dsharp"
        files[base]["solvertime"] = find_dsharp_time_cnt(f)
        files[base]["solverver"] = ["dsharp", "dsharp"]
        files[base]["confls"] = 0
        files[base]["comps"] = 0
    if ".out_minic2d" in f:
        files[base]["solver"] = "minic2d"
        files[base]["solvertime"] = find_minic2d_time_cnt(f)
        files[base]["solverver"] = ["minic2d", "minic2d"]
        files[base]["confls"] = 0
        files[base]["comps"] = 0
    if ".out_exactmc" in f:
        files[base]["solver"] = "exactmc"
        files[base]["solvertime"] = find_exactmc_time_cnt(f)
        files[base]["solverver"] = ["exactmc","exactmc"]
        files[base]["confls"] = 0
        files[base]["comps"] = 0


with open("mydata.csv", "w") as out:
    cols = "dirname,fname,"
    cols += "ganak_time,ganak_tout_t,ganak_mem_MB,ganak_call,ganak_ver,confls,comps"
    out.write(cols+"\n")
    for _, f in files.items():
        toprint = ""
        toprint += f["dirname"] + ","
        toprint += f["fname"] + ","

        # ganak_t
        if "solver" not in f:
            print("oops")
            print(f)
            continue

        # ganak_time
        if f["solvertime"] == [None, None]:
            toprint += ","
        else:
            toprint += "%s," % f["solvertime"][0]

        #ganak_tout_t, ganak_mem_MB, ganak_call
        if f["solvertout"] == [None, None]:
            toprint += ",,"
        else:
            toprint += "%s," % f["solvertout"][0]
            toprint += "%s," % f["solvertout"][1]
            toprint += "%s," % f["solvertout"][2]

        #ganak_ver
        if f["solverver"] == [None, None]:
            toprint += ",,"
        else:
            toprint += "%s-%s," % (f["solverver"][0], f["solverver"][1])
            toprint += "%s," % f["confls"]
            toprint += "%s"  % f["comps"]

        out.write(toprint+"\n")

os.system("sqlite3 mydb.sql < ganak.sqlite")
