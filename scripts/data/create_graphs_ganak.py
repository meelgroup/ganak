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
    res = cur.execute("SELECT ganak_ver FROM data where ganak_ver is not NULL and ganak_ver != '' group by ganak_ver")
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
todo = [ # "17e237e32f302198a52529900a4c58ee84a95ee0", # thursday night 22:07
        # "e10e390959c41990407439e16cec3b4ea8bfcf3a", # thu 17:1
        # "f6789ccc62c9748f03198467d2d24d2136901b1b", # thu 12:15
        # "5f3f40e7",
        # "7b5b5c09",
        # "89a6317b",
        # "db40b291",
        # "772a782c",
        "cbefa43a",
        # "cd4fa81d",
        # "a7031a8b",
        # "4b62829e", # this is really 73ab7a4992fce5e695355d6ac3d350cefac8d486
        # "3d72f8f8"
        ]
todo = versions
for ver in todo :
    dirs_call = get_dirs(ver)
    for dir,call in dirs_call:
        print("dir:", dir)
        print("call:", call)
        if "actexp 1.0" in call:
            continue
        # if "tdwithredbins 0" in call:
        #     continue
        # if "restart 1" in call:
        #     continue
        # if "vivif" not in call:
        #     continue
        if "probe 1" in call:
            continue
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
        fname2_s.append([fname2, call, ver[:6], num_solved, dir])

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
    f.write("plot [0:3600][:]\\\n")
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


os.system("gnuplot "+gnuplotfn)
print("okular run.eps")
os.system("okular run.eps")
