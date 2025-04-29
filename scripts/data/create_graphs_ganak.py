#!/bin/python3

import os
import sqlite3
import re
import sys

if len(sys.argv) < 2:
    print("ERROR: must call with --proj/--unproj/--all/--ganak/--example/--numbers")
    exit(-1)
fname_like = ""
# 20 examples
fname_like2 = " and fname in ('mc2023_track3_107.cnf', 'mc2024_track4_081.cnf', 'mc2023_track3_157.cnf', 'mc2023_track3_051.cnf', 'mc2023_track4_193.cnf', 'mc2023_track4_111.cnf', 'mc2024_track4_080.cnf', 'mc2023_track4_067.cnf', 'mc2024_track4_129.cnf', 'mc2023_track4_106.cnf', 'mc2023_track3_063.cnf', 'mc2023_track3_080.cnf', 'mc2023_track3_061.cnf', 'mc2023_track4_182.cnf', 'mc2024_track3_022.cnf', 'mc2024_track3_147.cnf', 'mc2023_track4_174.cnf', 'mc2023_track3_128.cnf', 'mc2024_track4_077.cnf', 'mc2024_track3_001.cnf', 'mc2023_track1_107.cnf', 'mc2024_track2-random_081.cnf', 'mc2023_track1_157.cnf', 'mc2023_track1_051.cnf', 'mc2023_track2_193.cnf', 'mc2023_track2_111.cnf', 'mc2024_track2-random_080.cnf', 'mc2023_track2_067.cnf', 'mc2024_track2-random_129.cnf', 'mc2023_track2_106.cnf', 'mc2023_track1_063.cnf', 'mc2023_track1_080.cnf', 'mc2023_track1_061.cnf', 'mc2023_track2_182.cnf', 'mc2024_track1_022.cnf', 'mc2024_track1_147.cnf', 'mc2023_track2_174.cnf', 'mc2023_track1_128.cnf', 'mc2024_track2-random_077.cnf', 'mc2024_track1_001.cnf', 'mc2023_track3_149.cnf') "
# 40 examples
# fname_like2 = "and fname in ('mc2023_track3_107.cnf', 'mc2024_track4_081.cnf', 'mc2023_track3_157.cnf', 'mc2023_track3_051.cnf', 'mc2023_track4_193.cnf', 'mc2023_track4_111.cnf', 'mc2024_track4_080.cnf', 'mc2023_track4_067.cnf', 'mc2024_track4_129.cnf', 'mc2023_track4_106.cnf', 'mc2023_track3_063.cnf', 'mc2023_track3_080.cnf', 'mc2023_track3_061.cnf', 'mc2023_track4_182.cnf', 'mc2024_track3_022.cnf', 'mc2024_track3_147.cnf', 'mc2023_track4_174.cnf', 'mc2023_track3_128.cnf', 'mc2024_track4_077.cnf', 'mc2024_track3_001.cnf', 'mc2024_track4_193.cnf', 'mc2024_track4_186.cnf', 'mc2023_track3_018.cnf', 'mc2023_track4_113.cnf', 'mc2024_track3_199.cnf', 'mc2023_track3_074.cnf', 'mc2024_track4_151.cnf', 'mc2023_track3_012.cnf', 'mc2024_track4_146.cnf', 'mc2024_track4_133.cnf', 'mc2023_track1_107.cnf', 'mc2024_track2-random_081.cnf', 'mc2023_track1_157.cnf', 'mc2023_track1_051.cnf', 'mc2023_track2_193.cnf', 'mc2023_track2_111.cnf', 'mc2024_track2-random_080.cnf', 'mc2023_track2_067.cnf', 'mc2024_track2-random_129.cnf', 'mc2023_track2_106.cnf', 'mc2023_track1_063.cnf', 'mc2023_track1_080.cnf', 'mc2023_track1_061.cnf', 'mc2023_track2_182.cnf', 'mc2024_track1_022.cnf', 'mc2024_track1_147.cnf', 'mc2023_track2_174.cnf', 'mc2023_track1_128.cnf', 'mc2024_track2-random_077.cnf', 'mc2024_track1_001.cnf', 'mc2024_track2-random_193.cnf', 'mc2024_track2-random_186.cnf', 'mc2023_track1_018.cnf', 'mc2023_track2_113.cnf', 'mc2024_track1_199.cnf', 'mc2023_track1_074.cnf', 'mc2024_track2-random_151.cnf', 'mc2023_track1_012.cnf', 'mc2024_track2-random_146.cnf', 'mc2024_track2-random_133.cnf', '2023_track3_149.cnf')"
if len(sys.argv) > 2 and sys.argv[2] == "--test":
  fname_like = fname_like2

if sys.argv[1] == "--example":
  with open("gen_table.sqlite", "w") as f:
    f.write(".mode table\n")
    only_dirs = [ "out-ganak-baseline",
                  "out-ganak-basic-sat-and-chronobt",
                  "out-ganak-also-enhanced-sat",
                  "out-ganak-also-dual-indep",
                  "out-ganak-also-extend-d-set" ]
    dirs = ""
    for dir in only_dirs:
      dirs += "'" + dir + "',"
    dirs = dirs[:-1]
    fname = "mc2023_track3_149.cnf"
    extra = ""
    f.write("select \
        replace(dirname,'out-ganak-mc','') as dirname,\
        (CASE WHEN ganak_time is not null then ganak_time else 9999999 END) as 'time', \
        ganak_mem_MB as 'mem(MB)', \
        (CASE WHEN ganak_time is NOT NULL THEN conflicts else NULL END)/(1000.0) as 'confls(K)', \
        (CASE WHEN ganak_time is not NULL THEN compsK else NULL END) as 'comps(K)',\
        (CASE WHEN ganak_time is not NULL THEN indep_sz else NULL END) as 'S-set',\
        (CASE WHEN opt_indep_sz!=indep_sz THEN opt_indep_sz else NULL END) as 'D-set'\
        from data where dirname IN ("+dirs+") and fname ='"+fname+"' order by time desc")
  os.system("sqlite3 mydb.sql < gen_table.sqlite")
  exit(0)

if sys.argv[1] == "--numbers":
    os.system("grep 'Init' out-gpmc/* > out-gpmc-init.txt")
    vars = []
    pvars = []
    with open("out-gpmc-init.txt", "r") as f:
      for line in f:
        line = line.strip()
        line = line.split()
        v = int(line[4])
        p = int(line[6][1:])
        pvars.append(p)
        vars.append(v)
    vars = sorted(vars)
    pvars = sorted(pvars)
    print("median vars: ", vars[len(vars)//2])
    print("median projected vars: ", pvars[len(pvars)//2])
    os.unlink("out-gpmc-init.txt")

    with open("gen_table.sqlite", "w") as f:
      f.write(".mode table\n")
      dir="out-ganak-also-extend-d-set"
      f.write("select 'data'");
      for col,col2 in [("indep_sz", "med S-set"), ("opt_indep_sz", "med D-set"), ("new_nvars", "med num vars after simp")]:
        f.write(", (SELECT "+col+" as 'median_"+col+"'\
        FROM data\
        where dirname IN ('"+dir+"') and "+col+" is not null"+fname_like+"\
        ORDER BY "+col+"\
        LIMIT 1\
        OFFSET (SELECT COUNT("+col+") FROM data\
          where dirname IN ('"+dir+"') \
          and "+col+" is not null) / 2) as '"+col2+"' \
      ")
      for col,col2 in [("gates_extended", "med syntactic ext."), ("padoa_extended", "med semantic ext.")]:
        f.write(", (SELECT "+col+" as 'median_"+col+"_NOZERO'\
        FROM data\
        where dirname IN ('"+dir+"') and "+col+" is not null "+fname_like+"\
                    and "+col+">0\
        ORDER BY "+col+"\
        LIMIT 1\
        OFFSET (SELECT COUNT("+col+") FROM data\
          where dirname IN ('"+dir+"') \
          and "+col+" is not null "+fname_like+" and "+col+">0) / 2) as '"+col2+"' \
      ")
      for col,col1, col2 in [("gates_extend_t", "gates_extended", "med syntactic ext t(s)"), ("padoa_extend_t", "padoa_extended", "med semantic ext t(s)")]:
        f.write(", (SELECT ROUND(avg("+col+"), 2) FROM data where dirname IN ('"+dir+"')  and "+col+" is not null) as '"+col2+"' ")
    os.system("sqlite3 mydb.sql < gen_table.sqlite")
    exit(0)


if sys.argv[1] != "--all" and sys.argv[1]!= "--proj" and sys.argv[1] != "--unproj" and sys.argv[1] != "--ganak":
    print("ERROR: must call with --proj/--unproj/--all/--ganak")
    exit(-1)

def convert_to_cactus(fname, fname2):
    # print("fname:" , fname)
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
not_calls = []
not_versions = []
only_calls = []
not_calls = []
todo = versions
if sys.argv[1] == "--unproj":
    fname_like = " and (fname like '%track1%' or fname like '%track2%') "
    only_dirs = [ "out-ganak-baseline",
                  "out-ganak-also-extend-d-set",
                  "out-gpmc",
                  "out-d4",
                  "out-sharptd"]
elif sys.argv[1] == "--proj":
    fname_like = " and (fname like '%track3%' or fname like '%track4%') "
    only_dirs = [ "out-ganak-baseline",
                  "out-ganak-also-extend-d-set",
                  "out-gpmc",
                  "out-d4"]
elif sys.argv[1] == "--all":
    fname_like = " "
    only_dirs = [ "out-ganak-baseline",
                  "out-ganak-also-extend-d-set",
                  "out-gpmc",
                  "out-d4"]
elif sys.argv[1] == "--ganak":
    fname_like = " "
    only_dirs = [ "out-ganak-baseline",
                  "out-ganak-basic-sat-and-chronobt",
                  "out-ganak-also-enhanced-sat",
                  "out-ganak-also-dual-indep",
                  "out-ganak-also-extend-d-set" ]
else:
    print("ERROR: must call with --proj/--unproj/--all")
    exit(-1)
if len(sys.argv) > 2 and sys.argv[2] == "--test":
  fname_like = fname_like + fname_like2

table_todo = []
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

        fname = "run-"+dir+".csv"
        with open("gencsv.sqlite", "w") as f:
            f.write(".headers off\n")
            f.write(".mode csv\n");
            f.write(".output "+fname+"\n")
            extra = ""
            f.write("select ganak_time from data where dirname='"+dir+"' and ganak_time is not NULL "+fname_like)
        os.system("sqlite3 mydb.sql < gencsv.sqlite")
        os.unlink("gencsv.sqlite")

        fname2 = fname + ".gnuplotdata"
        num_solved = convert_to_cactus(fname, fname2)
        fname2_s.append([fname2, call, ver[:10], num_solved, dir])
        table_todo.append([dir, ver])

if sys.argv[1] != "--ganak":
  with open("gen_table.sqlite", "w") as f:
    f.write(".mode table\n")
    # f.write(".mode colum\n")
    # f.write(".headers off\n")
    dirs = ""
    for dir,ver in table_todo:
      dirs += "'" + dir + "',"
    dirs = dirs[:-1]
    extra = ""
    f.write("select \
        replace(dirname,'out-ganak-mc','') as dirname,\
        sum(ganak_time is not null) as 'counted',\
        CAST(ROUND(sum(coalesce(ganak_time, 3600))/COUNT(*),0) AS INTEGER) as 'PAR2',\
        CAST(ROUND(avg(CASE WHEN ganak_time IS NOT NULL THEN ganak_mem_MB ELSE NULL END), 0) AS INTEGER) as 'avg mem(MB)', \
        sum(fname is not null) as 'nfiles'\
        from data where dirname IN ("+dirs+") "+fname_like+" group by dirname order by PAR2 desc")
  os.system("sqlite3 mydb.sql < gen_table.sqlite")
else:
  with open("gen_table.sqlite", "w") as f:
    f.write(".mode table\n")
    # f.write(".mode colum\n")
    # f.write(".headers off\n")
    dirs = ""
    for dir,ver in table_todo:
      dirs += "'" + dir + "',"
    dirs = dirs[:-1]
    extra = ""
    f.write("select \
        replace(dirname,'out-ganak-mc','') as dirname,\
        sum(ganak_time is not null) as 'counted',\
        CAST(ROUND(sum(coalesce(ganak_time, 3600))/COUNT(*),0) AS INTEGER) as 'PAR2',\
        CAST(ROUND(avg(CASE WHEN ganak_time IS NOT NULL THEN ganak_mem_MB ELSE NULL END), 0) AS INTEGER) as 'avg mem(MB)', \
        ROUND(avg(CASE WHEN ganak_time is NOT NULL THEN conflicts else NULL END)/(1000.0*1000.0), 2) as 'avg confls(M)', \
        ROUND(avg(CASE WHEN ganak_time is not NULL THEN compsK else NULL END)/(1000.0),2) as 'avg comps(M)',\
        sum(fname is not null) as 'nfiles'\
        from data where dirname IN ("+dirs+") "+fname_like+" group by dirname order by PAR2 desc")
  os.system("sqlite3 mydb.sql < gen_table.sqlite")

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
    # f.write("plot [:][10:]\\\n")
    f.write("plot [:][:]\\\n")
    i = 0
    towrite = ""
    for fn,call,ver,num_solved,dir in fname2_s:
        # if "restart" not in call and num_solved > 142:
        if True:
            call = re.sub("\"", "", call)
            dir  = re.sub("\"", "", dir)
            oneline = "\""+fn+"\" u 2:1 with linespoints  title \""+dir+"\""
            towrite += oneline
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
print("now please run: okular run.eps")
print("we'll run it now.. should pop up on your screen")
os.system("okular run.eps")
