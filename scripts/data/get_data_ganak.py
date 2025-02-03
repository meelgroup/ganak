#!/usr/bin/python3

import os
import glob
import decimal
import sys
import string

sys.set_int_max_str_digits(2000000)

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
            if "c [COUNTER] Elapsed time:" in line:
                t = float(line.split()[4])
            if "c exact arb int" in line:
                cnt = decimal.Decimal(line.split()[5]).log10()
            if "c s exact quadruple int" in line:
                cnt = decimal.Decimal(line.split()[5]).log10()

    # print("t:", t, "cnt: ", cnt)
    return [t,cnt]

#c o Components            = 339421
#c o conflicts             = 291050      (count 287680, sat 3370)
#c o decisions             = 643616      (count 335433, sat 308183)
def find_gpmc_time_cnt(fname):
    t = None
    cnt = None
    compsK = None
    conflicts = None
    decisionsK = None
    with open(fname, "r") as f:
        for line in f:
            line = line.strip()
            if "c o Real time " in line:
                t = float(line.split()[5])
            if "c s exact arb int" in line:
                if len(line.split()[5]) > 1000: cnt = len(line.split()[5])
                else: cnt = decimal.Decimal(line.split()[5]).log10()
            if "c s exact arb prec-sci" in line:
                if len(line.split()[5]) > 1000: cnt = len(line.split()[5])
                else: cnt = decimal.Decimal(line.split()[5]).log10()
            if "c o Components" in line:
              compsK = int(line.split()[4])/1000
            if "c o conflicts" in line:
              conflicts = int(line.split()[4])
            if "c o decisions" in line:
              decisionsK = int(line.split()[4])/1000

    return t,cnt,compsK,conflicts,decisionsK

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
            if "c s exact arb float" in line:
                cnt = decimal.Decimal(line.split()[5]).log10()

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
# c o Sampling set size: 94
# c o Opt sampling set size: 94
# c o CNF projection set size: xx
def find_arjun_time(fname):
    t = None
    backb_t = None
    backw_t = None
    unkn_sz = None
    orig_proj_sz = None
    indep_sz = None
    opt_indep_sz = None
    nvars = None
    gates_extend_t = None
    gates_extended = None
    padoa_extend_t = None
    padoa_extended = None
    with open(fname, "r") as f:
        for line in f:
            line = line.strip()
            if "c o Sampling set size:" in line:
              indep_sz = line.split()[5]
              indep_sz = indep_sz.strip()
              if indep_sz == "":
                indep_sz = 0
              else:
                indep_sz = int(indep_sz)
              if indep_sz == 4294967295:
                indep_sz = 0
            #c o opt ind size: 0 ind size: 0 nvars: 3463
            if "c o opt ind size" in line:
              nvars = int(line.split()[10])
            if "c o Opt sampling set size:" in line:
              opt_indep_sz = line.split()[6]
              opt_indep_sz = opt_indep_sz.strip()
              if opt_indep_sz == "":
                opt_indep_sz = 0
              else:
                opt_indep_sz = int(opt_indep_sz)
              if opt_indep_sz == 4294967295:
                opt_indep_sz = 0
            if "c o CNF projection set size:" in line:
              orig_proj_sz = int(line.split()[6])
            # c o [extend-gates] Gates added to opt indep: 26 T: 0.15
            # c o [arjun-extend] Start unknown size: 1027
            # c o [arjun-extend] Extend finished  orig size: 48 final size: 1056 Undef: 48 T: 25.36
            if "[extend-gates] Gates added to opt" in line:
              line = line.replace("[0m", "")
              line = line.replace("\x1b","")
              gates_extended = int(line.split()[8])
              gates_extend_t = float(line.split()[10])
              # print ("gates:", gates_extended, gates_extend_t)
            if "[arjun-extend] Extend finished" in line:
              line = line.replace("[0m", "")
              line = line.replace("\x1b","")
              orig = int(line.split()[7])
              final = int(line.split()[10])
              padoa_extended = final - orig
              padoa_extend_t = float(line.split()[14])
              # print ("padoa:", padoa_extended, padoa_extend_t)
            if "Start unknown size" in line:
              line = ''.join(filter(lambda x:x in string.printable, line))
              line = line.replace("[0m", "")
              unkn_sz = int(line.split()[6])
            if  "backward round finished" in line:
              line = ''.join(filter(lambda x:x in string.printable, line))
              line = line.replace("[0m", "")
              backw_t = float(line.split()[10])
            if  "% total" in line:
              if backb_t is None:
                backb_t = float(line.split()[2])
              else:
                backb_t += float(line.split()[2])
            if "c o Arjun T:" in line:
              assert t is None
              t = float(line.split()[4])
    return t, backb_t, backw_t, indep_sz, opt_indep_sz, orig_proj_sz, unkn_sz, nvars, gates_extended, gates_extend_t, padoa_extended, padoa_extend_t


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
            if "c o sat call/sat/unsat/confl/rst" in line:
              n = int(line.split()[4])
              rst = int(line.split()[8])
    return n,rst

#c o buddy called                   3
def find_bdd_called(fname):
    n = None
    with open(fname, "r") as f:
        for line in f:
            line = line.strip()
            if "c o buddy called /unsat ratio" in line:
              n = int(line.split()[6])
            elif "c o buddy called" in line:
              n = int(line.split()[4])
    return n

#c o Num restarts: 27
#c o [rst-cube] Num restarts: 1 orig cubes this rst: 0 total orig cubes: 0 total final cubes: 0 counted this rst: 0 total cnt so far: 0
def find_restarts(fname):
    n = None
    cubes_orig = 0
    cubes_final = 0
    with open(fname, "r") as f:
        for line in f:
            line = line.strip()

            # c o cubes orig/symm/final          1     / 0     / 1
            if "c o cubes orig:" in line:
              cubes_orig += int(line.split()[4])
              cubes_final += int(line.split()[6])
            if "c o Num restarts:" in line:
              n = int(line.split()[4])

            if "c o [rst-cube] Num restarts:" in line:
              cubes = int(line.split()[14])
              n = int(line.split()[5])
    return n,cubes_orig,cubes_final

#c o deletion done. T: 3.067
#c o cache pollutions call/removed  56828/28682
#c o cache K (lookup/ stores/ hits/ dels) 51     32     19     0       -- Klookup/s:  13.81
#c o cache pollutions call/removed  56828/28682
#c o cache miss rate                0.622
# c o [td] iter 189 best bag: 33 stepsK remain: -2 elapsed: 31.944
# c o [td] Primal graph   nodes: 548 edges: 8571 density: 0.029 edge/var: 15.669
def collect_cache_data(fname):
    cache_del_time = 0.0
    cache_miss_rate = None
    cache_lookupK = None
    cache_storeK = None
    density = None
    edge_var_ratio = None
    td_w = None
    td_t = None
    with open(fname, "r") as f:
        for line in f:
            line = line.strip()
            if "c o [td] iter" in line:
                td_w = int(line.split()[7])-1
                td_w = "%d" % td_w
                td_t = float(line.split()[12])
                td_t = "%f" % t
            if "c o [td] Primal graph" in line:
              density = float(line.split()[10])
              edge_var_ratio = float(line.split()[12])
            #c o cache miss rate                0.622
            if "c o cache miss rate" in line:
              cache_miss_rate = float(line.split()[5])
            #c o cache K (lookup/ stores/ hits/ dels) 51     32     19     0       -- Klookup/s:  13.81
            if "c o cache K (lookup/ stores/ hits/ dels)" in line:
              cache_lookupK = float(line.split()[8])
              cache_storeK = float(line.split()[9])
            if "c o deletion done. T:" in line:
              cache_del_time += float(line.split()[5])
    return cache_del_time, cache_miss_rate, cache_lookupK, cache_storeK, density, edge_var_ratio,td_w, td_t

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
                elif "doalarm 60" in call:
                  call = call.split("doalarm 60")[1]
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
    decisionsK = None
    with open(fname, "r") as f:
        for line in f:
            line = line.strip()
            if "c o decisions K" in line:
                decisionsK = int(line.split()[4])

    return decisionsK


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


def find_mem_out(fname):
    mem_out = 0
    with open(fname, "r") as f:
        for line in f:
            line = line.strip()
            if "std::bad_alloc" in line:
                mem_out = 1

    return mem_out


file_list = glob.glob("out-ganak-*/*cnf*")
file_list.extend(glob.glob("out-others-*/*cnf*"))
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

    if f.endswith(".timeout") or ".timeout_" in f:
        files[base]["solvertout"] = timeout_parse(f)

    files[base]["mem_out"] = find_mem_out(f)
    if  f.endswith(".out_ganak") or f.endswith(".out"):
        files[base]["solver"] = "ganak"
        files[base]["solvertime"] = find_ganak_time_cnt(f)
        files[base]["solverver"] = ganak_version(f)
        files[base]["conflicts"] = ganak_conflicts(f)
        files[base]["decisionsK"] = ganak_decisions(f)
        arjun_t, backb_t, backw_t, indep_sz, opt_indep_sz, orig_proj_sz, unkn_sz, new_nvars, gates_extended, gates_extend_t, padoa_extended, padoa_extend_t= find_arjun_time(f)
        files[base]["arjuntime"] = arjun_t
        files[base]["backbtime"] = backb_t
        files[base]["backwtime"] = backw_t
        files[base]["indepsz"] = indep_sz
        files[base]["optindepsz"] = opt_indep_sz
        files[base]["origprojsz"] = orig_proj_sz
        files[base]["unknsz"] = unkn_sz
        files[base]["newnvars"] = new_nvars
        files[base]["gates_extended"] = gates_extended
        files[base]["gates_extend_t"] = gates_extend_t
        files[base]["padoa_extended"] = padoa_extended
        files[base]["padoa_extend_t"] = padoa_extend_t
        cache_del_time, cache_miss_rate, cache_lookupK, cache_storeK, density, edge_var_ratio,td_w, td_t= collect_cache_data(f)
        files[base]["compsK"] = cache_lookupK
        files[base]["cache_miss_rate"] = cache_miss_rate
        files[base]["cache_storeK"] = cache_storeK
        files[base]["cache_del_time"] = cache_del_time
        files[base]["bddcalled"] = find_bdd_called(f)
        rst,cubes_orig,cubes_final = find_restarts(f)
        files[base]["restarts"] = rst
        files[base]["cubes_orig"] = cubes_orig
        files[base]["cubes_finale"] = cubes_final
        sat_called,sat_rst = find_sat_called(f)
        files[base]["sat_called"] = sat_called
        files[base]["satrst"] = sat_rst
        files[base]["td_width"] = td_w
        files[base]["td_time"] = td_t
        files[base]["primal_density"] = density
        files[base]["primal_edge_var_ratio"] = edge_var_ratio
    if ".out_approxmc" in f:
        files[base]["solver"] = "approxmc"
        files[base]["solvertime"] = find_approxmc_time_cnt(f)
        files[base]["solverver"] = approxmc_version(f)
    if ".out_gpmc" in f:
        files[base]["solver"] = "gpmc"
        t,cnt,compsK,conflicts,decisionsK = find_gpmc_time_cnt(f)
        files[base]["solvertime"] = [t, cnt]
        files[base]["conflicts"] = conflicts
        files[base]["compsK"] = compsK
        files[base]["decisionsK"] = decisionsK
        files[base]["solverver"] = ["gpmc", "gpmc"]
    if ".out_sharptd" in f:
        files[base]["solver"] = "sharptd"
        files[base]["solvertime"] = find_sharpsat_time_cnt(f)
        files[base]["solverver"] = ["sharptd", "sharptd"]
        td = sstd_treewidth(f)
        files[base]["td_width"] = td[0]
        files[base]["td_time"] = td[1]
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
    cols += "ganak_time,ganak_tout_t,ganak_mem_MB,ganak_call,ganak_ver,conflicts,decisionsK,compsK,primal_density,primal_edge_var_ratio,td_width,td_time,arjun_time,backboneT,backwardT,indepsz,optindepsz,origprojsz,new_nvars,unknsz,cache_del_time,cache_miss_rate,bdd_called,sat_called,sat_rst,rst,cubes_orig,cubes_final,mem_out,gates_extended,gates_extend_t,padoa_extended,padoa_extend_t"
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

        if "conflicts" not in f or f["conflicts"] is None:
            toprint += ","
        else:
          toprint += "%s," % f["conflicts"]

        if "decisionsK" not in f or f["decisionsK"] is None:
            toprint += ","
        else:
          toprint += "%d," % f["decisionsK"]

        if "compsK" not in f or f["compsK"] is None:
            toprint += ","
        else:
          toprint += "%s," % f["compsK"]

        if "primal_density" not in f or f["primal_density"] is None:
            toprint += ","
        else:
          toprint += "%s," % f["primal_density"]

        if "primal_edge_var_ratio" not in f or f["primal_edge_var_ratio"] is None:
            toprint += ","
        else:
          toprint += "%s," % f["primal_edge_var_ratio"]

        if "td_width" not in f or f["td_width"] is None:
            toprint += ","
        else:
          toprint += "%s," % f["td_width"]

        if "td_time" not in f or f["td_time"] is None:
            toprint += ","
        else:
          toprint += "%s,"  % f["td_time"]

        if "arjuntime" not in f or f["arjuntime"] is None:
            toprint += ","
        else:
          toprint += "%s,"  % f["arjuntime"]

        if "backbtime" not in f or f["backbtime"] is None:
            toprint += ","
        else:
          toprint += "%s,"  % f["backbtime"]

        if "backwtime" not in f or f["backwtime"] is None:
            toprint += ","
        else:
          toprint += "%s,"  % f["backwtime"]

        if "indepsz" not in f or f["indepsz"] is None:
            toprint += ","
        else:
          toprint += "%s,"  % f["indepsz"]

        if "optindepsz" not in f or f["optindepsz"] is None:
            toprint += ","
        else:
          toprint += "%s,"  % f["optindepsz"]

        if "origprojsz" not in f or f["origprojsz"] is None:
            toprint += ","
        else:
          toprint += "%s,"  % f["origprojsz"]

        if "newnvars" not in f or f["newnvars"] is None:
            toprint += ","
        else:
          toprint += "%s,"  % f["newnvars"]

        if "unknsz" not in f or f["unknsz"] is None:
            toprint += ","
        else:
          toprint += "%s,"  % f["unknsz"]

        if "cache_del_time" not in f or f["cache_del_time"] is None:
            toprint += ","
        else:
          toprint += "%s,"  % f["cache_del_time"]

        if "cache_miss_rate" not in f or f["cache_miss_rate"] is None:
            toprint += ","
        else:
          toprint += "%s,"  % f["cache_miss_rate"]

        if "bddcalled" not in f or f["bddcalled"] is None:
            toprint += ","
        else:
          toprint += "%s,"  % f["bddcalled"]

        if "sat_called" not in f or f["sat_called"] is None:
            toprint += ","
        else:
          toprint += "%s,"  % f["sat_called"]

        if "satrst" not in f or f["satrst"] is None:
            toprint += ","
        else:
          toprint += "%s,"  % f["satrst"]

        if "restarts" not in f or f["restarts"] is None:
            toprint += ","
        else:
          toprint += "%s,"  % f["restarts"]

        if "cubes_orig" not in f or f["cubes_orig"] is None:
            toprint += ","
        else:
          toprint += "%s,"  % f["cubes_orig"]

        if "cubes_final" not in f or f["cubes_final"] is None:
            toprint += ","
        else:
          toprint += "%s,"  % f["cubes_final"]

        if "mem_out" not in f or f["mem_out"] is None:
            toprint += ","
        else:
          toprint += "%s,"  % f["mem_out"]

        if "gates_extended" not in f or f["gates_extended"] is None:
            toprint += ","
        else:
          toprint += "%s,"  % f["gates_extended"]

        if "gates_extend_t" not in f or f["gates_extend_t"] is None:
            toprint += ","
        else:
          toprint += "%s,"  % f["gates_extend_t"]

        if "padoa_extended" not in f or f["padoa_extended"] is None:
            toprint += ","
        else:
          toprint += "%s,"  % f["padoa_extended"]

        if "padoa_extend_t" not in f or f["padoa_extend_t"] is None:
            toprint += ""
        else:
          toprint += "%s"  % f["padoa_extend_t"]

        out.write(toprint+"\n")

os.system("sqlite3 mydb.sql < ganak.sqlite")
