#!/bin/python3

import argparse
import os
import sqlite3
import re
import nbformat as nbf
import plotext as plt


def convert_to_cactus(fname, fname2):
    with open(fname, "r") as f:
        times = sorted(float(line.split()[0]) for line in f)

    with open(fname2, "w") as f2:
        for i, t in enumerate(times):
            f2.write(f"{i+1} \t{t}\n")
    return len(times)


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


def build_csv_data(todo, only_dirs, only_calls, not_calls, not_versions, fname_like, verbose=False):
    fname2_s = []
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
                if verbose:
                    print(f"  Skipping dir={dir} ver={ver}")
                continue

            if verbose:
                print(f"  Processing dir={dir} ver={ver}")
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
    return fname2_s, table_todo


def print_summary_tables(table_todo, fname_like, full=False, verbose=False):
    dirs = ",".join("'" + dir + "'" for dir, _ in table_todo)
    vers = ",".join("'" + ver + "'" for _, ver in table_todo)

    compact_cols = [
        ("replace(dirname,'out-ganak-mc','')",                       "dirname"),
        ("replace(ganak_call,'././ganak_','')",                      "call"),
        ("sum(ganak_time is not null)",                              "solved"),
        ("CAST(ROUND(sum(coalesce(ganak_time,3600))/COUNT(*),0) AS INTEGER)", "PAR2"),
        ("ROUND(avg(conflicts)/(1000.0*1000.0), 2)",                 "av confM"),
        ("CAST(ROUND(avg(ganak_mem_mb),0) AS INTEGER)",              "av memMB"),
        ("sum(signal == 11)",                                        "sigSEGV"),
        ("sum(signal == 6)",                                         "sigABRT"),
        ("sum(signal == 14)",                                        "sigALRM"),
        ("sum(signal == 8)",                                         "sigFPE"),
        ("ROUND(avg(cache_miss_rate),2)",                            "av cmiss"),
    ]
    full_only_cols = [
        ("sum(mem_out)",                                             "mem out"),
        ("ROUND(avg(decisionsK)/(1000.0), 2)",                       "av decM"),
        ("CAST(ROUND(avg(sat_called/1000.0),0) AS INTEGER)",         "av satcK"),
        ("CAST(avg(opt_indep_sz-indep_sz) AS INTEGER)",              "avg-diff-opt-sz"),
        ("ROUND(avg(gates_extend_t), 3)",                            "gates-ext-t"),
        ("ROUND(avg(padoa_extend_t), 3)",                            "padoa-ext-t"),
        ("ROUND(avg(gates_extended), 3)",                            "gates-ext"),
        ("ROUND(avg(padoa_extended), 3)",                            "padoa-ext"),
        ("CAST(ROUND(max(cache_del_time), 0) AS INTEGER)",           "max cachdT"),
        ("CAST(ROUND(avg(backbone_time),0) AS INTEGER)",             "av backT"),
        ("CAST(ROUND(avg(arjun_time),0) AS INTEGER)",                "av arjT"),
        ("CAST(ROUND(avg(td_time),0) AS INTEGER)",                   "av tdT"),
        ("ROUND(avg(td_width),0)",                                   "av tdw"),
        ("ROUND(avg(cache_avg_hit_vars),2)",                         "av chitvs"),
        ("ROUND(avg(compsK/1000.0),2)",                              "av compsM"),
        ("sum(fname is not null)",                                   "nfiles"),
    ]

    cols = compact_cols + (full_only_cols if full else [])
    select_clause = ",\n        ".join(f"{expr} as '{alias}'" for expr, alias in cols)

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
            query = (f"select\n        {select_clause}\n"
                     f"        from data where dirname IN ({dirs}) and ganak_ver IN ({vers})"
                     f" {fname_like} {counted_req}group by dirname order by solved asc")
            if verbose:
                print(f"  Summary query: {query[:120]}...")
            f.write(query)
        os.system("sqlite3 mydb.sql < gen_table.sqlite")
        os.unlink("gen_table.sqlite")


def _median_subquery(col, dir, ver, fname_like, nozero=False):
    extra = f" and {col}>0" if nozero else ""
    base = f"dirname='{dir}' and ganak_ver='{ver}' and {col} is not null{fname_like}{extra}"
    return (f"(SELECT {col} FROM data WHERE {base}"
            f" ORDER BY {col} LIMIT 1"
            f" OFFSET (SELECT COUNT({col}) FROM data WHERE {base}) / 2)")


def print_median_tables(table_todo, fname_like, verbose=False):
    if not table_todo:
        return

    plain_cols  = ["indep_sz", "opt_indep_sz", "orig_proj_sz", "new_nvars", "ganak_mem_mb"]
    nozero_cols = ["gates_extended", "padoa_extended"]

    union_parts = []
    for dir, ver in table_todo:
        parts = [f"'{dir}'", f"'{ver}'"]
        for col in plain_cols:
            parts.append(f"{_median_subquery(col, dir, ver, fname_like)} as median_{col}")
        for col in nozero_cols:
            parts.append(f"{_median_subquery(col, dir, ver, fname_like, nozero=True)} as median_{col}_NOZERO")
        union_parts.append("select " + ", ".join(parts))

    query = "\nUNION ALL\n".join(union_parts)
    if verbose:
        print(f"  Median table query ({len(table_todo)} rows)")
    with open("gen_table.sqlite", "w") as f:
        f.write(".mode table\n")
        f.write(query + "\n")
    os.system("sqlite3 mydb.sql < gen_table.sqlite")
    os.unlink("gen_table.sqlite")


def print_distribution(table_todo, fname_like, col, label, xscale="linear", xmin=None, xlabel=None):
    for dir, ver in table_todo:
        con = sqlite3.connect("mydb.sql")
        cur = con.cursor()
        res = cur.execute(
            "SELECT " + col + " FROM data WHERE dirname='" + dir +
            "' AND ganak_ver='" + ver +
            "' AND " + col + " IS NOT NULL"
            " AND (ganak_time - arjun_time) >= 100" + fname_like
        )
        values = [row[0] for row in res]
        con.close()

        if not values:
            print(f"No {label} data for {dir}")
            continue

        title = f"{label}: {dir} [ganak_time-arjun_time >= 100s, n={len(values)}]"
        print(f"\n=== {title} ===")
        plt.clf()
        plt.theme("dark")
        plt.plot_size(160, 30)
        if xmin is not None:
            values = [v for v in values if v >= xmin]
        if xscale == "log":
            import math
            values = [math.log10(v) for v in values if v > 0]
        plt.hist(values, bins=20)
        if xmin is not None:
            plt.xlim(xmin, max(values))
        plt.xlabel(xlabel if xlabel is not None else col)
        plt.ylabel("count")
        plt.show()


def print_sigabrt_files(table_todo, fname_like):
    dirs = ",".join("'" + dir + "'" for dir, _ in table_todo)
    vers = ",".join("'" + ver + "'" for _, ver in table_todo)
    con = sqlite3.connect("mydb.sql")
    cur = con.cursor()
    cur.execute(
        f"SELECT COUNT(*) FROM data WHERE dirname IN ({dirs}) AND ganak_ver IN ({vers})"
        f" AND signal=6{fname_like}"
    )
    count = cur.fetchone()[0]
    if count == 0:
        con.close()
        return
    print(f"\n::: WARNING: {count} instance(s) with sigABRT (signal=6) :::")
    cur.execute(
        f"SELECT dirname, fname, timeout_t FROM data WHERE dirname IN ({dirs}) AND ganak_ver IN ({vers})"
        f" AND signal=6{fname_like} ORDER BY dirname, fname"
    )
    rows = cur.fetchall()
    con.close()
    str_rows = [(d, f, f"{t:.2f}" if t is not None else "N/A") for d, f, t in rows]
    widths = [max(len(h), max(len(r[i]) for r in str_rows)) for i, h in enumerate(("dirname", "fname", "timeout_t"))]
    sep = "+-" + "-+-".join("-" * w for w in widths) + "-+"
    fmt = "| " + " | ".join(f"{{:<{w}}}" for w in widths) + " |"
    print(sep)
    print(fmt.format("dirname", "fname", "timeout_t"))
    print(sep)
    for row in str_rows:
        print(fmt.format(*row))
    print(sep)


def print_distributions(table_todo, fname_like):
    print_distribution(table_todo, fname_like, "cache_miss_rate",  "cache miss rate")
    print_distribution(table_todo, fname_like, "compsK",           "num components (K) [log10 x-axis]", xscale="log", xmin=1, xlabel="LOG compsK")
    print_distribution(table_todo, fname_like, "td_width",         "TD width", xmin=0)
    print_distribution(table_todo, fname_like, "ganak_mem_mb",     "memory usage (MB) [log10 x-axis]", xscale="log", xmin=1, xlabel="LOG mem_mb")


def generate_gnuplot(fname2_s, verbose=False):
    gnuplotfn = "run-all.gnuplot"
    if verbose:
        print(f"Writing gnuplot script to {gnuplotfn} with {len(fname2_s)} data series")
    with open(gnuplotfn, "w") as f:
        f.write("set term postscript eps color lw 1 \"Helvetica\" 12 size 9,5\n")
        f.write("set output \"run.eps\"\n")
        f.write("set title \"Counter ganak\"\n")
        f.write("set notitle\n")
        f.write("set key bottom right\n")
        # f.write("set xtics 200\n")
        f.write("set logscale x\n")
        f.write("unset logscale y\n")
        f.write("set ylabel  \"Instances counted\"\n")
        f.write("set xlabel \"Time (s)\"\n")
        # f.write("plot [:][10:]\\\n")
        # f.write("plot [500:4000][1000:1200]\\\n")
        f.write("plot [0.1:3600][0.1:]\\\n")
        # f.write(" \"runkcbox-prearjun.csv.gnuplotdata\" u 2:1 with linespoints  title \"KCBox\",\\\n")
        # f.write(" \"runsharptd-prearjun.csv.gnuplotdata\" u 2:1 with linespoints  title \"SharptTD\",\\\n")
        towrite = ""
        for fn, call, ver, _, dir in fname2_s:
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
    return gnuplotfn


def create_notebook(dirs):
    nb = nbf.v4.new_notebook()
    texts = []

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

    cells = [nbf.v4.new_code_cell(t) for t in texts]
    nb['cells'] = cells

    filename = 'overview.ipynb'
    with open(filename, 'w') as f:
        nbf.write(nb, f)
    print(f"Notebook '{filename}' created successfully.")


# ---- Configuration ----

versions = get_versions()
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
    # "out-ganak-mccomp2324-1152658-1", # fixing bug, testing pura and oraclemult
    # "out-ganak-mccomp2324-1193808-0", # fixing counting bug, adding new cube extend system, fixing cube counting (and maybe effectiveness)
    "out-ganak-mccomp2324-1229753-0", # lots of bug fixes, beauty changes with Claude, etc
    "out-ganak-mccomp2324-1231407-0", # the same as above but without (most) of the Claude improvements
]
# only_dirs = [ "mei-march-2026-1239767" ]

# not_calls = ["--nvarscutoffcache 20", "--nvarscutoffcache 3"]
# not_calls = ["--satsolver 0"]
not_versions = []
# only_dirs = []
# only_calls = ["vsadsadjust 128"]
only_calls = []
# not_calls = ["restart"]
not_calls = []
todo = versions

# fname filter examples (pass via --fname on the command line):
# unweighted only: --fname "%track1%" "%track3%"
# weighted only:   --fname "%track2%" "%track4%"
# unproj only:     --fname "%track1%" "%track2%"
# proj only:       --fname "%track3%" "%track4%"


# ---- Main ----

def main():
    parser = argparse.ArgumentParser(description="Generate cactus plots and tables for ganak benchmark data")
    parser.add_argument("--verbose", "-v", action="store_true", help="Print progress information")
    parser.add_argument("--full", action="store_true", help="Print full summary table (default: compact)")
    parser.add_argument("--nograph", action="store_true", help="Skip opening the graph in okular")
    parser.add_argument("--fname", nargs="+", metavar="PATTERN", default=[],
                        help="Filter by fname pattern(s), e.g. --fname '%%track1%%' '%%track3%%'")
    args = parser.parse_args()

    if args.fname:
        clauses = " or ".join(f"fname like '{p}'" for p in args.fname)
        fname_like = f" and ({clauses}) "
    else:
        fname_like = ""

    if args.verbose:
        print(f"Found {len(versions)} versions in database")
        print("Building CSV data...")
    fname2_s, table_todo = build_csv_data(todo, only_dirs, only_calls, not_calls, not_versions, fname_like, args.verbose)

    if args.verbose:
        print(f"Selected {len(table_todo)} dir/version combinations")
    seen = set()
    for dir, _ in table_todo:
        if dir not in seen:
            seen.add(dir)
            os.system(f"./cache_miss_bucket_summary.py {dir}")

    print_distributions(table_todo, fname_like)

    if args.verbose:
        print("Printing summary tables...")
    print_summary_tables(table_todo, fname_like, args.full, args.verbose)
    print_sigabrt_files(table_todo, fname_like)

    if args.verbose:
        print("Printing median tables...")
    print_median_tables(table_todo, fname_like, args.verbose)

    if args.verbose:
        print("Generating gnuplot script...")
    gnuplotfn = generate_gnuplot(fname2_s, args.verbose)

    dirs = ",".join("'" + dir + "'" for dir, _ in table_todo)
    create_notebook(dirs)

    for path in ["run.eps", "run.pdf", "run.png"]:
        if os.path.exists(path):
            os.unlink(path)
    os.system("gnuplot "+gnuplotfn)
    os.system("epstopdf run.eps run.pdf")
    os.system("pdftoppm -png run.pdf run")
    if not args.nograph:
        print("okular run.eps")
        os.system("okular run.eps")


if __name__ == "__main__":
    main()
