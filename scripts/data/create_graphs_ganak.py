#!/bin/python3

import os
import sqlite3
import re
import nbformat as nbf


def convert_to_cactus(fname, fname2):
    with open(fname, "r") as f:
        time = [float(line.split()[0]) for line in f]

    lastnum = -1
    with open(fname2, "w") as f2:
        for a in range(0, 3600, 1):
            num = sum(1 for t in time if t < a)
            if lastnum != num:
                f2.write(f"{num} \t{a}\n")
            lastnum = num
    return len(time)


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


def get_dirs(ver: str):
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


def gnuplot_name_cleanup(name: str) -> str:
    # remove all non-alphanumeric characters except for underscores and dashes
    name = re.sub(r'\"', '', name)
    # replace multiple underscores or dashes with a single one
    name = re.sub(r'_', '=', name)
    return name


versions = get_versions()
fname2_s = []
# not_calls = ["ExactMC"]
# exactm: out-ganak-6318929.pbs101-5/
# exactmc2: out-ganak-6328707.pbs101-7
# sharpsat: out-ganak-6318929.pbs101-7

only_dirs = [
    # "out-ganak-mccomp2324-14299825",
    # "out-ganak-mccomp2324-14309534-0", # this IS the current best
    # "out-ganak-mccomp2324-14309534-1", # messed up
    # "out-ganak-mccomp2324-14321462-" # let's try blocking without messing it up
    # "out-ganak-mccomp2324-14362931-0", # new hash, timestamp archetype
    # "out-ganak-mccomp2324-14391524-", # changeable recompute archetype cutoff, better cache usage
    # "out-ganak-mccomp2324-14513393-0", # update arjun's order and add gate-based-eqlit, also fix gate-based-rem by not doing strengthening
    # "out-ganak-mccomp2324-14513746-0", # back to previous, except the strengthening fix for gate-based-rem
    # "out-ganak-mccomp2324-14519572-", # as previous 2, but fixed parsing of stupid weights
    # "out-ganak-mccomp2324-14528572-", # binaries in holder are no longer updated
    # "out-ganak-mccomp2324-14538562-1", # full-probe before everything in remove_eq_lits
    # "out-ganak-mccomp2324-14549970-1", # LTO
    # "out-ganak-mccomp2324-14552736-", # no LTO, but some change (e.g. sandybridge etc)
    # "out-ganak-mccomp2324-14558670-1", # some stuff enabled previously for no flto
    # "out-ganak-mccomp2324-14575065-", # occ-bve-resolv after bve
    # "out-ganak-mccomp2324-14576195-", # as 2 above, but now playing with TD iters
    # "out-ganak-mccomp2324-14621165", # float, also try different TD exp
    # "out-ganak-mccomp2324-14621994-0", # float fixed mem free, try different tdmaxw
    # "out-ganak-mccomp2324-14624577", # exact, no cache check, exp check
    # "out-ganak-mccomp2324-14624580", # exact, more exp check
    # "out-ganak-mccomp2324-14637456", # try vivif options
    # "out-ganak-mccomp2324-14648650-0", # much more precise cache
    # "out-ganak-mccomp2324-14655977-", # get bins from cadiback
    # "out-ganak-mccomp2324-14656347-0", # get bins, scc, enable weakening, more weakening -- GOOD
    # "out-ganak-mccomp2324-14661402-", # try different weakenings
    # "out-ganak-mccomp2324-14675861-0", # submitted to mccomp 2025 -- GOOD
    # "out-ganak-mccomp2324-14675861-", # default setup along WITHOUT appmc. Trying tditers -- ganak_7d97636055e_9104724fa_26d64aac --tdlookiters 20
    # "out-ganak-mccomp2324-14675861-1", # default setup along WITH appmc. Trying tditers
    # "out-ganak-mccomp2324-15010600-0", # TD start from 0
    # "out-ganak-mccomp2324-3382-0", # TRILLIUM -- TD start from 0
    # "out-ganak-mccomp2324-21238", # TRILLIUM --non-eq, and eq, and no probabilistic
    # "out-ganak-mccomp2324-21349-0", # TRILLIUM, ganak_7d97636055e_9104724fa_26d64aac (i.e. old run that was the fastest)
    # "out-ganak-mccomp2324-21481-0", # TRILLIUM, fixing memory usage for --prob 0
    # "out-ganak-mccomp2324-33393-0", # TRILLIUM, fixing td starting from 0, and maybe other TD issues too
    # "out-ganak-mccomp2324-35541-0", # TRILLIUM, fixing td starting from 0, now cutting disjoint components at toplevel for correct centroid
    # "out-ganak-mccomp2324-37704-0", # TRILLIUM, trying less scc upfront, weird non-committed code
    # "out-ganak-mccomp2324-41122-0", # TRILLIUM, new TD score, but BAD, not updating depth when all vars have been decided
    # "out-ganak-mccomp2324-41923-0", # TRILLIUM, new TD score, fixed depth update finally
    # "out-ganak-mccomp2324-43391-0", # TRILLIUM, new TD score, fixed depth update
    # "out-ganak-mccomp2324-44137-0", # TRILLIUM, finally good fix for TD
    # "out-ganak-mccomp2324-983755", # TRILLIUM, new arjun with oracle cache fix, autarky
    # "out-ganak-mccomp2324-984574", # same as above, bug not sigFPE and no loop in autarky
    # "out-ganak-mccomp2324-1140000-0", # dodgy cache update
    # "out-ganak-mccomp2324-1140184-", # no dodgy cache, check SBVA
    # "out-ganak-mccomp2324-1146702-", # sbva checks, oracle mult checks
    "out-ganak-mccomp2324-1152658-1", # fixing bug, testing pura and oraclemult
    "out-ganak-mccomp2324-1193808-0", # fixing counting bug, adding new cube extend system, fixing cube counting (and maybe effectiveness)
    "out-ganak-mccomp2324-1229753-0", # lots of bug fixes, beauty changes with Claude, etc
    "out-ganak-mccomp2324-1231407-0", # the same as above but without (most) of the Claude improvements
]
# only_dirs = [ "mei-march-2026-1237508-" ]

# not_calls = ["--nvarscutoffcache 20", "--nvarscutoffcache 3"]
# not_calls = ["--satsolver 0"]
not_versions = []
# only_dirs = []
# only_calls = ["vsadsadjust 128"]
only_calls = []
# not_calls = ["restart"]
not_calls = []
todo = versions

# all
fname_like = ""
# unweighted only
# fname_like = " and (fname like '%track1%' or fname like '%track3%') "
# weighted only
# fname_like = " and (fname like '%track2%' or fname like '%track4%') "
# unproj only
# fname_like = " and (fname like '%track1%' or fname like '%track2%') "
# proj only
# fname_like = " and (fname like '%track3%' or fname like '%track4%') "

table_todo = []
for ver in todo:
    dirs_call = get_dirs(ver)
    for dir, call in dirs_call:
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
            if not inside:
                bad = True

        if len(only_dirs) != 0:
            inside = False
            for only_dir in only_dirs:
                if only_dir in (dir+"/"):
                    inside = True
            if not inside:
                bad = True

        if bad:
            continue
        fname = "run-"+dir+".csv"
        with open("gencsv.sqlite", "w") as f:
            f.write(".headers off\n")
            f.write(".mode csv\n")
            f.write(".output "+fname+"\n")
            f.write("select ganak_time from data where dirname='"+dir+"' and ganak_ver='"+ver+"'\n and ganak_time is not NULL "+fname_like)
        os.system("sqlite3 mydb.sql < gencsv.sqlite")
        os.unlink("gencsv.sqlite")

        fname2 = fname + ".gnuplotdata"
        num_solved = convert_to_cactus(fname, fname2)
        fname2_s.append([fname2, call, ver[:10], num_solved, dir])
        table_todo.append([dir, ver])

for only_counted in [False, True]:
    counted_req = ""
    if only_counted:
        print("::: --------- Data based on ONLY benchmarks that are COUNTED ------- :::")
        counted_req = " and ganak_time is not NULL "
    else:
        print("::: --------- Data based on ALSO UNCOUNTED benchmarks ------- :::")
    with open("gen_table.sqlite", "w") as f:
        f.write(".mode table\n")
        # f.write(".mode colum\n")
        # f.write(".headers off\n")
        dirs = ""
        vers = ""
        for dir, ver in table_todo:
            dirs += "'" + dir + "',"
            vers += "'" + ver + "',"
        dirs = dirs[:-1]
        vers = vers[:-1]
        f.write("select \
        replace(dirname,'out-ganak-mc','') as dirname,\
        replace(ganak_call,'././ganak_','') as call,\
        sum(mem_out) as 'mem out', \
        sum(signal == 11) as 'sigSEGV', \
        sum(signal == 6) as 'sigABRT', \
        sum(signal == 14) as 'sigALRM', \
        sum(signal == 8) as 'sigFPE', \
        CAST(ROUND(avg(ganak_mem_MB), 0) AS INTEGER) as 'av memMB',\
        ROUND(avg(conflicts)/(1000.0*1000.0), 2) as 'av confM', \
        ROUND(avg(decisionsK)/(1000.0), 2) as 'av decM', \
        CAST(ROUND(avg(sat_called/1000.0),0) AS INTEGER) as 'av satcK',\
        sum(ganak_time is not null) as 'solved',\
        CAST(ROUND(sum(coalesce(ganak_time, 3600))/COUNT(*),0) AS INTEGER) as 'PAR2',\
        CAST(avg(opt_indep_sz-indep_sz) AS INTEGER) as 'avg-diff-opt-sz',\
        ROUND(avg(gates_extend_t), 3) as 'gates-ext-t',\
        ROUND(avg(padoa_extend_t), 3) as 'padoa-ext-t',\
        ROUND(avg(gates_extended), 3) as 'gates-ext',\
        ROUND(avg(padoa_extended), 3) as 'padoa-ext',\
        CAST(ROUND(max(cache_del_time), 0) AS INTEGER) as 'max cachdT',\
        CAST(ROUND(avg(backbone_time),0) AS INTEGER) as 'av backT',\
        CAST(ROUND(avg(arjun_time),0) AS INTEGER) as 'av arjT',\
        CAST(ROUND(avg(td_time),0) AS INTEGER) as 'av tdT',\
        ROUND(avg(td_width),0) as 'av tdw',\
        ROUND(avg(cache_miss_rate),2) as 'av cmiss',\
        ROUND(avg(cache_avg_hit_vars),2) as 'av chitvs',\
        ROUND(avg(ganak_mem_mb),2) as 'av memMB',\
        ROUND(avg(compsK/1000.0),2) as 'av compsM',\
        sum(fname is not null) as 'nfiles'\
        from data where dirname IN ("+dirs+") and ganak_ver IN ("+vers+") "+fname_like+" "+counted_req+"group by dirname order by solved asc")
    os.system("sqlite3 mydb.sql < gen_table.sqlite")
    os.unlink("gen_table.sqlite")

for dir, ver in table_todo:
    with open("gen_table.sqlite", "w") as f:
        f.write(".mode table\n")
        f.write("select '"+dir+"', '"+ver+"'")
        for col in "indep_sz", "opt_indep_sz", "orig_proj_sz", "new_nvars", "ganak_mem_mb":
            f.write(", (SELECT "+col+" as 'median_"+col+"'\
        FROM data\
        where dirname IN ('"+dir+"') and ganak_ver IN ('"+ver+"') and "+col+" is not null"+fname_like+"\
        ORDER BY "+col+"\
        LIMIT 1\
        OFFSET (SELECT COUNT("+col+") FROM data\
          where dirname IN ('"+dir+"') and ganak_ver IN ('"+ver+"') \
          and "+col+" is not null) / 2) as median_"+col+" \
      ")
        for col in "gates_extended", "padoa_extended":
            f.write(", (SELECT "+col+" as 'median_"+col+"_NOZERO'\
        FROM data\
        where dirname IN ('"+dir+"') and ganak_ver IN ('"+ver+"') and "+col+" is not null "+fname_like+"\
                    and "+col+">0\
        ORDER BY "+col+"\
        LIMIT 1\
        OFFSET (SELECT COUNT("+col+") FROM data\
          where dirname IN ('"+dir+"') and ganak_ver IN ('"+ver+"') \
          and "+col+" is not null "+fname_like+" and "+col+">0) / 2) as median_"+col+" \
      ")
    os.system("sqlite3 mydb.sql < gen_table.sqlite")
    os.unlink("gen_table.sqlite")

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
    # f.write("plot [500:4000][1000:1200]\\\n")
    f.write("plot [:][:]\\\n")
    # f.write(" \"runkcbox-prearjun.csv.gnuplotdata\" u 2:1 with linespoints  title \"KCBox\",\\\n")
    # f.write(" \"runsharptd-prearjun.csv.gnuplotdata\" u 2:1 with linespoints  title \"SharptTD\",\\\n")
    towrite = ""
    for fn, call, ver, num_solved, dir in fname2_s:
        # if "restart" not in call and num_solved > 142:
        if True:
            call = gnuplot_name_cleanup(call)
            dir  = gnuplot_name_cleanup(dir)
            ver  = gnuplot_name_cleanup(ver)
            oneline = "\""+fn+"\" u 2:1 with linespoints  title \""+ver+"-"+dir+"-"+call+"\""
            towrite += oneline
            towrite += ",\\\n"
    towrite = towrite[:(len(towrite)-4)]
    f.write(towrite)

for path in ["run.eps", "run.pdf", "run.png"]:
    if os.path.exists(path):
        os.unlink(path)


def create_notebook():
    # Create a new notebook
    nb = nbf.v4.new_notebook()
    texts = []

    # Define the text content for the markdown cells
    text = """
# Step 1: Import necessary libraries
import pandas as pd
import sqlite3
import matplotlib.pyplot as plt
from functools import reduce
import numpy as np

dirs = ["""+dirs+"""]
"""
    texts.append(text)

    colnames = ["ganak_time", "ganak_mem_MB", "conflicts", "decisionsK",
        "compsK", "td_width", "arjun_time", "backbone_time", "indep_sz", "td_time", "sat_called", "cache_miss_rate"]
    for colname in colnames:
        text = """
colname='"""+colname+"""'
names=[]
dfs = []
conn = sqlite3.connect('mydb.sql')
for d in dirs:
  # Step 3: Run the SQL query and load the results into a DataFrame
  query = "SELECT fname, "+colname+", ganak_call FROM data where "+colname+" is not NULL and dirname='"+d+"' order by "+colname
  df1 = pd.read_sql_query(query, conn)
  df1['num'] = range(len(df1))
  dfs.append(df1)
  names.append(d+" " +df1['ganak_call'][0])

for i in range(len(dfs)):
    for c in dfs[i].columns:
      if c == 'num': continue
      dfs[i].rename(columns={c: c+'_'+str(i)}, inplace=True)
result = dfs[0]

for i, df in enumerate(dfs[1:], start=1):
    result = pd.merge(result, df, on='num', how='outer', suffixes=(f'_{i-1}', f'_{i}'))

# Merge all DataFrames
conn.close()

# Step 5: Plot the data
plt.figure(figsize=(10, 6))
for i in range(len(dirs)):
  values = result[(colname + '_' + str(i))]
  min_non_zero = values[values > 0].min()
  adjusted_values = np.where(values == 0, min_non_zero, values)
  log_values = np.log10(adjusted_values)
  plt.plot(result['num'],log_values,marker='o')
plt.title('Plot of CNFs counted vs. log10 '+colname)
plt.xlabel('Number of CNFs counted')
plt.ylabel('log10 '+colname)
plt.legend(names,loc='center left', bbox_to_anchor=(0, -0.3))
plt.grid(True)
plt.show()
  """
        texts.append(text)

    for col1, col2 in [("td_width", "ganak_time"), ("td_width", "td_time"), ("td_width", "arjun_time"), ("td_width", "backbone_time"), ("compsK", "ganak_time")]:
        text = """
col1='"""+col1+"""'
col2='"""+col2+"""'
dirs = ["""+dirs+"""]

names = []
dfs = []
conn = sqlite3.connect('mydb.sql')

# Assign a color for each dirname
colors = plt.cm.tab10(np.linspace(0, 1, len(dirs)))

for d in dirs:
    # Run the SQL query and load the results into a DataFrame
    query = f"SELECT fname, {col1}, {col2}, ganak_call FROM data WHERE {col1} IS NOT NULL AND {col2} is not NULL and dirname='{d}' ORDER BY {colname}"
    df1 = pd.read_sql_query(query, conn)
    dfs.append(df1)
    names.append(d + " " + df1['ganak_call'][0])

conn.close()

# Plot the data
plt.figure(figsize=(10, 10))
for i, d in enumerate(dirs):
    plt.scatter(np.log10(dfs[i][col1]), np.log10(dfs[i][col2]), color=colors[i], label=names[i])

plt.title('Scatterplot')
plt.xlabel(f"{col1} (log10 scale)")
plt.ylabel(f"{col2} (log10 scale)")
plt.legend(loc='center left', bbox_to_anchor=(0, -0.2))
plt.show()
"""
        texts.append(text)

    # Create code cells
    cells = [nbf.v4.new_code_cell(t) for t in texts]
    nb['cells'] = cells

    filename = 'overview.ipynb'
    with open(filename, 'w') as f:
        nbf.write(nb, f)
    print(f"Notebook '{filename}' created successfully.")


create_notebook()
os.system("gnuplot "+gnuplotfn)
os.system("epstopdf run.eps run.pdf")
os.system("pdftoppm -png run.pdf run")
print("okular run.eps")
os.system("okular run.eps")
