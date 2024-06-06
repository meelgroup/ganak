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
            if "c exact arb int" in line:
                cnt = decimal.Decimal(line.split()[5]).log10()

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
            if "s pmc" in line:
                cnt = decimal.Decimal(line.split()[2])

    return [t,cnt]

#c o Arjun T: 206.14
def find_arjun_time(fname):
    t = None
    backb_t = None
    with open(fname, "r") as f:
        for line in f:
            line = line.strip()
            if line[:3] == "c o" and "% total " in line:
              if backb_t is None:
                backb_t = float(line.split()[2])
              else:
                backb_t+ = float(line.split()[2])
            if "c o Arjun T:" in line:
              assert t is None
              t = float(line.split()[4])
    return t, backb_t

#c o sat call/sat/unsat/conflK/rst  0     0     0     0     0
#c o sat called/sat/unsat/conflK    6     6     0     0
def find_sat_called(fname):
    n = None
    rst = None
    with open(fname, "r") as f:
        for line in f:
            line = line.strip()
            if "c o sat called/sat/unsat/conflK" in line:
              n = int(line.split()[4])
            if "c o sat call/sat/unsat/conflK/rst" in line:
              n = int(line.split()[4])
              rst = int(line.split()[8])
    return n,rst

#c o buddy called                   3
def find_bdd_called(fname):
    n = None
    with open(fname, "r") as f:
        for line in f:
            line = line.strip()
            if "c o buddy called" in line:
              n = int(line.split()[4])
    return n

#c o Num restarts: 27
#c o [rst-cube] Num restarts: 1 orig cubes this rst: 0 total orig cubes: 0 total final cubes: 0 counted this rst: 0 total cnt so far: 0
def find_restarts(fname):
    n = None
    cubes = 0
    with open(fname, "r") as f:
        for line in f:
            line = line.strip()

            if "c o Total num cubes:" in line:
              cubes += int(line.split()[5])
            if "c o Num restarts:" in line:
              n = int(line.split()[4])

            if "c o [rst-cube] Num restarts:" in line:
              cubes = int(line.split()[14])
              n = int(line.split()[5])
    return n,cubes

#c o deletion done. T: 3.067
def collect_cache_deletion_time(fname):
    t = 0.0
    with open(fname, "r") as f:
        for line in f:
            line = line.strip()
            if "c o deletion done. T:" in line:
              t += float(line.split()[5])
    return t

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
                  call = " ".join(call.split()[:-1])
                call = call.replace(" -t real", "")
                if "doalarm 3600" in call:
                  call = call.split("doalarm 3600")[1]
                else:
                  call = call.split("doalarm 900")[1]
                call = call.replace("././ganak ", "")
                call = call.replace("././d4-1d9cc6146f18b8 ", "")
                call = call.replace("././approxmc ", "")
                call = call.replace("././gpmc ", "")
                call = call.replace("././KCBox-371eb601f2aa", "")
                call = call.replace("././sharpSAT", "")
                call = call.strip()

    return [t, m, call]

# c o [td] iter 189 best bag: 33 stepsK remain: -2 elapsed: 31.944
def ganak_treewidth(fname) -> list[str]:
    tw = ""
    t = ""
    with open(fname, "r") as f:
        for line in f:
            line = line.strip()
            if "c o [td] iter" in line:
                tw = int(line.split()[7])-1
                tw = "%d" % tw
                t = float(line.split()[12])
                t = "%f" % t
    return [tw, t]

# c o width 45
# c o CMD: timeout 60.000000s ./flow_cutter_pace17 <tmp/instance1709332682043115_14040_59941_1.tmp >tmp/instance1709332682043118_14040_59941_2.tmp 2>/dev/null
def sstd_treewidth(fname) -> list[str]:
    tw = ""
    t = ""
    with open(fname, "r") as f:
        for line in f:
            line = line.strip()
            if "c o width" in line:
                tw = int(line.split()[3])
                tw = "%d" % tw
            if "CMD: timeout" in line:
                t = float(line.split()[4][:-1])
                t = "%f" % t
    return [tw, t]


# c o conflicts                      2503077      -- confl/s:    1434.45
def ganak_conflicts(fname) -> str:
    confl = ""
    with open(fname, "r") as f:
        for line in f:
            line = line.strip()
            if "c o conflicts" in line:
                confl = int(line.split()[3])
                confl = "%d" % confl

    return confl


#c o decisions K                    8         -- Kdec/s:     0.71
def ganak_decisions(fname) -> int|None:
    decisions = None
    with open(fname, "r") as f:
        for line in f:
            line = line.strip()
            if "c o decisions K" in line:
                decisions = int(line.split()[4])

    return decisions


# c o cache K (lookup/ stores/ hits) 67229  34594  32634   -- Klookup/s:  38.53
def ganak_comps(fname) -> str:
    comps = ""
    with open(fname, "r") as f:
        for line in f:
            line = line.strip()
            if "c o cache K (lookup/ stores/ hits)" in line:
                comps = int(line.split()[7])
                comps = "%d" % comps
            if "c o cache K (lookup/ stores/ hits/ dels)" in line:
                comps = int(line.split()[8])
                comps = "%d" % comps

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
            if "c o CMS revision" in line:
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
        files[base]["decisions"] = ganak_decisions(f)
        files[base]["comps"] = ganak_comps(f)
        arjun_t, backb_t = find_arjun_time(f)
        files[base]["arjuntime"] = arjun_t
        files[base]["backbtime"] = backb_t
        files[base]["cachedeltime"] = collect_cache_deletion_time(f)
        files[base]["bddcalled"] = find_bdd_called(f)
        rst,cubes = find_restarts(f)
        files[base]["restarts"] = rst
        files[base]["cubes"] = cubes
        sat_called,sat_rst = find_sat_called(f)
        files[base]["satcalled"] = sat_called
        files[base]["satrst"] = sat_rst
        td = ganak_treewidth(f)
        files[base]["td-width"] = td[0]
        files[base]["td-time"] = td[1]
    if ".out_approxmc" in f:
        files[base]["solver"] = "approxmc"
        files[base]["solvertime"] = find_approxmc_time_cnt(f)
        files[base]["solverver"] = approxmc_version(f)
    if ".out_gpmc" in f:
        files[base]["solver"] = "gpmc"
        files[base]["solvertime"] = find_gpmc_time_cnt(f)
        files[base]["solverver"] = ["gpmc", "gpmc"]
    if ".out_sharptd" in f:
        files[base]["solver"] = "sharptd"
        files[base]["solvertime"] = find_sharpsat_time_cnt(f)
        files[base]["solverver"] = ["sharptd", "sharptd"]
        td = sstd_treewidth(f)
        files[base]["td-width"] = td[0]
        files[base]["td-time"] = td[1]
    if ".out_d4" in f:
        files[base]["solver"] = "d4"
        files[base]["solvertime"] = find_d4_time_cnt(f)
        files[base]["solverver"] = ["d4", "d4"]
    if ".out_dsharp" in f:
        files[base]["solver"] = "dsharp"
        files[base]["solvertime"] = find_dsharp_time_cnt(f)
        files[base]["solverver"] = ["dsharp", "dsharp"]
    if ".out_minic2d" in f:
        files[base]["solver"] = "minic2d"
        files[base]["solvertime"] = find_minic2d_time_cnt(f)
        files[base]["solverver"] = ["minic2d", "minic2d"]
    if ".out_exactmc" in f:
        files[base]["solver"] = "exactmc"
        files[base]["solvertime"] = find_exactmc_time_cnt(f)
        files[base]["solverver"] = ["exactmc","exactmc"]


with open("mydata.csv", "w") as out:
    cols = "dirname,fname,"
    cols += "ganak_time,ganak_tout_t,ganak_mem_MB,ganak_call,ganak_ver,confls,decs,comps,td_width,td_time,arjun_time,backbone_time,cache_del_time,bdd_called,sat_called,sat_rst,rst,cubes"
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
        if f["solvertime"][0] is None:
            toprint += ","
        else:
            toprint += "%s," % f["solvertout"][0]
            # toprint += "%s," % f["solvertime"][0]

        #ganak_tout_t, ganak_mem_MB, ganak_call
        if f["solvertout"] == [None, None, None]:
            toprint += ",,,"
        else:
            toprint += "%s," % f["solvertout"][0]
            toprint += "%s," % f["solvertout"][1]
            toprint += "%s," % f["solvertout"][2]

        #ganak_ver
        if f["solverver"] == [None, None]:
            toprint += ","
        else:
          toprint += "%s-%s," % (f["solverver"][0], f["solverver"][1])

        if "confls" not in f or f["confls"] is None:
            toprint += ","
        else:
          toprint += "%s," % f["confls"]

        if "decisions" not in f or f["decisions"] is None:
            toprint += ","
        else:
          toprint += "%d," % f["decisions"]

        if "comps" not in f or f["comps"] is None:
            toprint += ","
        else:
          toprint += "%s," % f["comps"]

        if "td-width" not in f or f["td-width"] is None:
            toprint += ","
        else:
          toprint += "%s," % f["td-width"]

        if "td-time" not in f or f["td-time"] is None:
            toprint += ","
        else:
          toprint += "%s,"  % f["td-time"]

        if "arjuntime" not in f or f["arjuntime"] is None:
            toprint += ","
        else:
          toprint += "%s,"  % f["arjuntime"]

        if "backbtime" not in f or f["backbtime"] is None:
            toprint += ","
        else:
          toprint += "%s,"  % f["backbtime"]

        if "cachedeltime" not in f or f["cachedeltime"] is None:
            toprint += ","
        else:
          toprint += "%s,"  % f["cachedeltime"]

        if "bddcalled" not in f or f["bddcalled"] is None:
            toprint += ","
        else:
          toprint += "%s,"  % f["bddcalled"]

        if "satcalled" not in f or f["satcalled"] is None:
            toprint += ","
        else:
          toprint += "%s,"  % f["satcalled"]

        if "satrst" not in f or f["satrst"] is None:
            toprint += ","
        else:
          toprint += "%s,"  % f["satrst"]

        if "restarts" not in f or f["restarts"] is None:
            toprint += ","
        else:
          toprint += "%s,"  % f["restarts"]

        if "cubes" not in f or f["cubes"] is None:
            toprint += ""
        else:
          toprint += "%s"  % f["cubes"]


        out.write(toprint+"\n")

os.system("sqlite3 mydb.sql < ganak.sqlite")
