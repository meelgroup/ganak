#!/bin/python3

import os
import sqlite3
import re


def convert_to_cactus(fname, fname2):
    print("fname:" , fname)
    f2 = open(fname2, "w")
    f = open(fname, "r")
    text = f.read()
    mylines = text.splitlines()
    i = 0;
    time = []
    for line in mylines:
      time.append(float(line.split()[0]))
      i += 1

    lastnum = -1
    for a in range(0, 3600, 1):
      num = 0
      for t in time:
        #print "t: %f a: %d" %(t, a)
        if (t < a) :
          num += 1

      if (lastnum != num):
          f2.write("%d \t%d\n" %(num, a))
      lastnum = num
    f.close()
    f2.close()
    return len(mylines)


def get_versions():
    vers = []
    con = sqlite3.connect("mydb.sql")
    cur = con.cursor()
    res = cur.execute("""
                      SELECT ganak_ver
                      FROM data
                      where ganak_ver is not NULL and ganak_ver != '' group by ganak_ver""")
    for a in res:
        vers.append(a[0])

    con.close()
    return vers

def get_dirs(ver : str):
    ret = []
    con = sqlite3.connect("mydb.sql")
    cur = con.cursor()
    res = cur.execute("SELECT dirname, ganak_call FROM data where ganak_ver='"+ver+"' group by dirname")
    for a in res:
        call = a[1]
        call = re.sub("././ganak", "", call)
        call = re.sub(" mc2022.*cnf.*", "", call)
        ret.append([a[0], call])

    con.close()
    return ret

versions = get_versions()
fname2_s = []
# not_calls = ["ExactMC"]
not_calls = []
# not_versions = ["sharpsat", "gpmc", "6368237b"]
# only_calls = ["--ignore 1 --arjun 1 --maxcache 3500 --vivif 1 --decide 2 --sbva 1000"]
# not_versions = ["ganak"]
# not_calls = ["forcebranch", "target"] # "cachetime"
# exactm: out-ganak-6318929.pbs101-5/
# exactmc2: out-ganak-6328707.pbs101-7
# sharpsat: out-ganak-6318929.pbs101-7

only_dirs = [
  # "6683080",
             # "out-ganak-6683080.pbs101-1/", # best of GANAK so far, lbd 1
             #"out-ganak-6396805.pbs101-1/",  # best, but no cache
        #     "out-ganak-6396805.pbs101-12/", # used to be best
        #     "out-ganak-6318929.pbs101-4/", # approxmc
             #"out-ganak-6749880.pbs101-0/", # like out-ganak-6683080.pbs101-1, but cluster is slower(!!)
             # "out-ganak-6683080.pbs101-4/", # last run, without refactored Arjun+fixed SAT
             #"out-ganak-6841576.pbs101-0/", # trying TD exponent, turning off gates for Arjun
             #"out-ganak-6853697.pbs101-4/", # backbone first, if proj is smaller than TD, td=0.1, turn off gates -- ran with CMS 58286e58ef64156ed2c16ec4104820c7fb4aeda5
             #"out-ganak-6867174.pbs101-0", # fix with max grow 6, also try no freq
             #"out-ganak-6870225.pbs101-0", # max grow 6, maxtd 5
        #     "out-ganak-6896385.pbs101-0", # supposedly new, not really exciting
             # "out-ganak-6907885.pbs101-", # BAD variable activity merged idea
             # "out-ganak-6909164.pbs101-", # BAD var activity merged, polarity changed
             #"out-ganak-6916780.pbs101-", # variable activity
             #"out-ganak-6910211.pbs101-2", # check diferent component sortings
            # "out-ganak-6917417.pbs101-", # sortings, lbd, polarity
            # "out-ganak-6918163.pbs101-0", # restart on mccomp 2023, 20k restarts
             # "out-ganak-6965806.pbs101", # mccomp 2022+3, restarts base etc
             #"out-ganak-7015279.pbs101", # new extension and symm run
  #######
             # "out-ganak-7021521.pbs101-0", # 2023 mccomp, restarts, fixed blocking lit, cube extend only no contract
             #"out-ganak-7028560.pbs101-", # 2022 + 2023
             #"out-ganak-7041485.", # kr-24 instances from Arijit
             # "out-ganak-7041554.pbs101-", # rapid restarts
             # "out-ganak-7047972.pbs101-0",
             # "out-ganak-7049225.pbs101-",
             # "out-ganak-7022833.pbs101-4", # best ever
             # "out-ganak-7048280.pbs101-0", # best ever now

             # "out-ganak-7173534.pbs101-0", # proj-2023 first run, 4GB
             # "out-ganak-7180435.pbs101-0", # proj-2023 4B, vsads
             # "out-ganak-7184202.pbs101-" #proj 2023, 4GB, freq-128 (vsads), checking vsads params
             # "out-ganak-7189867.pbs101-1" # proj 2023, 4GB, cache elem simplification, TD back to 0.1, SAT restarts
             # "out-ganak7197939",
             # "out-ganak-7206369.pbs101-",
             # TODO: higher maxw than out-ganak-7205692.pbs101-6, rerun fixed cadiback
             # "out-ganak-7178163.pbs101-0", # proj-2023 16GB, ganak
             # "out-ganak-7178163.pbs101-2", # proj-2023 16GB, d4
             # "out-ganak-7178163.pbs101-3", # proj-2023 16GB, gpmc
             # "out-ganak-7205692.pbs101-6", # old best
             # #"out-ganak-7247003.pbs101-", # fixed cadiback, --tdmaxw 40 seems best, but 15 is fine too actually # BEST
             # "out-ganak-7255014.pbs101-0", # different tdminw-s, also try lbd 2. LBD is indifferent, tdminw high is BAD
             # "out-ganak-7266814.pbs101-1", # contract over TDW, higher tdmaxw, higher tdminw. Was buggy in a few ways. Let's re-run. Fixed contraction in the meanwhile.
             # "out-ganak-7294423.pbs101-0", # fixed memout from contraction, fixed too much extend, distill-bin, BIG before backbone
             # "out-ganak-7307327.pbs101-4", # fixed cadical, fixed resolv-subs, fixed oracle, 2 new options for extra oracle & resolv-subs
             #"out-ganak-7308235.pbs101-", # newest cadical -- ALL very good
                                             # best is "--tdmaxw 10 --varfreqdiv 30"
             # BAD PARAMS.... "out-ganak-7316071.pbs101-", # also appmc, play around with some options
             # "out-ganak-7320764.pbs101-0", #best: --tdmaxw 10 --varfreqdiv 25 , --appmct 2000
             # "out-ganak-7320764.pbs101-9", #best: --tdmaxw 10 --varfreqdiv 25 , --appmct 2000
             # "out-ganak-7320764.pbs101-0", #best: --tdmaxw 10 --varfreqdiv 25 , --appmct 2000 (could try varfreqdiv 20)
             "out-ganak-7334726.pbs101-1", #  new run, some parameter tuning, fixed SAT restart bug (again, this time good) -- BEST, non-appmc, --tdmaxw 10 --varfreqdiv 25
             "out-ganak-7334726.pbs101-8", #  new run, some parameter tuning, fixed SAT restart bug (again, this time good) -- BEST, appmc
             #"out-ganak-7348395.pbs101-" # better extend -- maybe extend is not good?
  # todo: no sbva, no extend, --compsort 5 & 4,
             "out-ganak-7366311.pbs101" # --backbone 0 is good, --tdexpmult 0.3 + --tdmaxw 100 is good


             # "out-ganak-6318929.pbs101-5/", # exactmc
             # "out-ganak-6328707.pbs101-7/", # exactmc
             # "out-ganak-6318929.pbs101-7/", # sharpsat
             # no point in combining out-ganak-7178422.pbs101-0 with out-ganak-7184237.pbs101-2, 171 either way
             # "out-ganak-7178422.pbs101-2", # unproj-2023 16 GB d4
             # "out-ganak-7178422.pbs101-3", # unproj-2023 16 GB gpmc
             # "out-ganak-7178422.pbs101-", # unproj-2023 16 GB
             # "out-ganak-7180395.pbs101-", #unproj 2023 4GB, freq-128 (vsads)
             # "out-ganak-7184237.pbs101-2", # unproj-2023 4GB, vsads, checking vsads params
             # "out-ganak-7189905.pbs101-" # unproj 2023, 4GB, cache elem simplification, TD back to 0.1, SAT restarts
             # "out-ganak-7206355.pbs101", # unproj, new freq, not terrible with high TW weight
             # "out-ganak-7225470.pbs101-6",
             # "out-ganak-7225470.pbs101-7",
             # "out-ganak-7162995.pbs101-0", # new run, good. -- old freq setup
             # "out-ganak-7246958", # fixed cadiback, high tdmaxw
             # "out-ganak-7255018.pbs101-7", # tdmaxw 90 + higher tdminw. 12 minw seems to work well -- BEST (maxw 90, minw 12)
             # "out-ganak-7257776.pbs101-9", # higher tdmaxw, higher tdminw -- BEST (maxw 100, minw 15)
             #"out-ganak-7295009.pbs101-0", #fixed memout from contraction, fixed too much extend, distill-bin, BIG before backbone
             # "out-ganak-7306993.pbs101-0",# fixed cadical, fixed resolv-subs, fixed oracle, 2 new options for extra oracle & resolv-subs
                                          # best is "--tdminw 15 --tdmaxw 100 --resolvsub 0 --extraoracle 1"
             # BAD PARAMS "out-ganak-7316065.pbs101", # also appmc, play around with some opitions
             # "out-ganak-7320968.pbs101-1", # best: --tdminw 15 --tdmaxw 100, --appmct 2000
             # "out-ganak-7320968.pbs101-7", # best: --tdminw 15 --tdmaxw 100, --appmct 2000
             # "out-ganak-7320968.pbs101-1", # best: --tdminw 15 --tdmaxw 100 --varfreqdiv 30, --appmct 2000 (tdexpmul 1?)
             # "out-ganak-7320968.pbs101-7", # best: --tdminw 15 --tdmaxw 100 --varfreqdiv 30, --appmct 2000 (tdexpmul 1?)
             # "out-ganak-7334807.pbs101-2", # new run, some parameter tuning, fixed SAT restart bug (again, this time good) BEST non-appmc
             # "out-ganak-7334807.pbs101-7", # new run, some parameter tuning, fixed SAT restart bug (again, this time good) BEST appmc


             ]
# only_dirs = ["out-ganak-6828273"] #-- functional synth
#"6393432", "6393432", "6349002",, "6349002", "6387743" "6356951"] #, "out-ganak-6318929.pbs101-4", "out-ganak-6328707.pbs101-7", "out-ganak-6318929.pbs101-7"] #,"6348728" "6346880", "6335522", "6328982", "6328707"]
# "6349002",
# only_dirs = ["6606250"]
# not_calls = ["--nvarscutoffcache 20", "--nvarscutoffcache 30", "--nvarscutoffcache 40", "--nvarscutoffcache 1", "--nvarscutoffcache 2",  "--nvarscutoffcache 3"]
not_versions = []
# only_calls = ["--lbd 1"] #
# only_dirs = []
# only_calls = ["--polar"]
only_calls = ["compsort"]
only_calls = []
todo = versions
for ver in todo :
    dirs_call = get_dirs(ver)
    for dir,call in dirs_call:
        bad = False
        for not_call in not_calls:
          if not_call in call:
            bad = True
        for not_version in not_versions:
          if not_version in ver:
            bad = True

        if len(only_calls) != 0:
          inside = False
          for only_call in only_calls:
            if only_call in call:
              inside = True
          if not inside: bad = True

        if len(only_dirs) != 0:
          inside = False
          for only_dir in only_dirs:
            if only_dir in (dir+"/"):
              inside = True
          if not inside: bad = True

        if bad:
          continue
        print("----")
        print("dir:", dir)
        print("call:", call)
        print("ver:", ver)

        # if "actexp 1.0" in call:
        #     continue
        # if "tdwithredbins 0" in call:
        #     continue
        # if "restart 1" in call:
        #     continue
        # if "vivif" not in call:
        #     continue
        # if "probe 1" in call:
        #     continue
        # if "polar" not in call:
        #     continue
        fname = "run-"+dir+".csv"
        with open("gencsv.sqlite", "w") as f:
            f.write(".headers off\n")
            f.write(".mode csv\n");
            f.write(".output "+fname+"\n")
            f.write("select ganak_time from data where dirname='"+dir+"' and ganak_ver='"+ver+"'\n and ganak_time is not NULL")
        os.system("sqlite3 mydb.sql < gencsv.sqlite")
        os.unlink("gencsv.sqlite")

        fname2 = fname + ".gnuplotdata"
        num_solved = convert_to_cactus(fname, fname2)
        fname2_s.append([fname2, call, ver[:10], num_solved, dir])

gnuplotfn = "run-all.gnuplot"
with open(gnuplotfn, "w") as f:
    f.write("set term postscript eps color lw 1 \"Helvetica\" 12 size 9,5\n")
    f.write("set output \"run.eps\"\n")
    f.write("set title \"Counter ganak\"\n")
    f.write("set notitle\n")
    f.write("set key bottom right\n")
    # f.write("set xtics 200\n")
    f.write("unset logscale x\n")
    f.write("unset logscale y\n")
    f.write("set ylabel  \"Instances counted\"\n")
    f.write("set xlabel \"Time (s)\"\n")
    f.write("plot [:][120:]\\\n")
    i = 0
    # f.write(" \"runkcbox-prearjun.csv.gnuplotdata\" u 2:1 with linespoints  title \"KCBox\",\\\n")
    # f.write(" \"runsharptd-prearjun.csv.gnuplotdata\" u 2:1 with linespoints  title \"SharptTD\",\\\n")
    towrite = ""
    for fn,call,ver,num_solved,dir in fname2_s:
        # if "restart" not in call and num_solved > 142:
        if True:
            towrite +="\""+fn+"\" u 2:1 with linespoints  title \""+ver+"-"+dir+"-"+call+"\""
            towrite +=",\\\n"
    towrite = towrite[:(len(towrite)-4)]
    f.write(towrite)


if os.path.exists("run.eps"):
  os.unlink("run.eps")
if os.path.exists("run.pdf"):
  os.unlink("run.pdf")
if os.path.exists("run.png"):
  os.unlink("run.png")

os.system("gnuplot "+gnuplotfn)
os.system("epstopdf run.eps run.pdf")
os.system("pdftoppm -png run.pdf run")
print("okular run.eps")
os.system("okular run.eps")
