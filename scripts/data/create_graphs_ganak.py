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
        #     "out-ganak-6318929.pbs101-5/", # exactmc
        #     "out-ganak-6328707.pbs101-7/", # exactmc
             #"out-ganak-6396805.pbs101-1/",  # best, but no cache
        #     "out-ganak-6396805.pbs101-12/", # used to be best
        #     "out-ganak-6318929.pbs101-4/", # approxmc
        #     "out-ganak-6318929.pbs101-7/", # sharpsat
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
             # "out-ganak-7041554.pbs101-0", # rapid restarts
             # "out-ganak-7047972.pbs101-0",
             # "out-ganak-7049225.pbs101-",
             #"out-ganak-7022833.pbs101-4", # best ever
             #"out-ganak-7048280.pbs101-0", # best ever now
             #"out-ganak-7162995.pbs101-", # new run, good.
             "out-ganak-7173534.pbs101", # proj-2023 first run
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
    f.write("plot [599:][:]\\\n")
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
