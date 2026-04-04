#!/bin/python3

import argparse
import base64
import itertools
import math
import os
import sqlite3
import re
import nbformat as nbf

BLUE   = "\033[94m"
RED    = "\033[91m"
RESET = "\033[0m"


def convert_to_cactus(fname, fname2):
    with open(fname, "r") as f:
        times = sorted(float(line.split()[0]) for line in f)

    with open(fname2, "w") as f2:
        for i, t in enumerate(times):
            f2.write(f"{i+1} \t{t}\n")
    return len(times)


def get_versions():
    vers = []
    con = sqlite3.connect("data.sqlite3")
    cur = con.cursor()
    res = cur.execute("""
                      SELECT ganak_ver
                      FROM data
                      where ganak_ver is not NULL and ganak_ver != '' group by ganak_ver""")
    for a in res:
        vers.append(a[0])
    con.close()
    return vers


def get_matching_dirs(only_dirs):
    """Return all dirnames from DB prefixed by any entry in only_dirs.
    Returns all dirnames if only_dirs is empty."""
    con = sqlite3.connect("data.sqlite3")
    cur = con.cursor()
    res = cur.execute("SELECT DISTINCT dirname FROM data")
    all_dirs = [row[0] for row in res]
    con.close()
    if not only_dirs:
        return all_dirs
    return [d for d in all_dirs if any((d + "/").startswith(p) for p in only_dirs)]


def get_dirs(ver: str):
    ret = []
    con = sqlite3.connect("data.sqlite3")
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


def build_csv_data(todo, matched_dirs, only_calls, not_calls, not_versions, fname_like, verbose=False):
    fname2_s = []
    table_todo = []
    matched_dirs_set = set(matched_dirs)
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

            if dir not in matched_dirs_set:
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
            os.system("sqlite3 data.sqlite3 < gencsv.sqlite")
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
            title = "Data based on ONLY benchmarks that are COUNTED"
        else:
            title = "Data based on ALSO UNCOUNTED benchmarks"
        print(f"\n{BLUE}{title}{RESET}")
        if only_counted:
            counted_req = " and ganak_time is not NULL "
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
        os.system("sqlite3 data.sqlite3 < gen_table.sqlite")
        os.unlink("gen_table.sqlite")


def _median_subquery(col, dir, ver, fname_like, nozero=False):
    extra = f" and {col}>0" if nozero else ""
    base = f"dirname='{dir}' and ganak_ver='{ver}' and {col} is not null{fname_like}{extra}"
    return (f"(SELECT {col} FROM data WHERE {base}"
            f" ORDER BY {col} LIMIT 1"
            f" OFFSET (SELECT COUNT({col}) FROM data WHERE {base}) / 2)")


def _avg_subquery(col, dir, ver, fname_like):
    base = f"dirname='{dir}' and ganak_ver='{ver}' and {col} is not null{fname_like}"
    return f"(SELECT CAST(ROUND(AVG({col}),0) AS INTEGER) FROM data WHERE {base})"


def print_median_tables(table_todo, fname_like, verbose=False):
    if not table_todo:
        return

    title = "Median values per directory"
    print(f"\n{BLUE}{title}{RESET}")

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
    os.system("sqlite3 data.sqlite3 < gen_table.sqlite")
    os.unlink("gen_table.sqlite")


def print_instance_stats_table(table_todo, fname_like, verbose=False):
    if not table_todo:
        return

    title = "Instance stats: variables, independent set sizes, irredundant clauses (median / avg)"
    print(f"\n{BLUE}{title}{RESET}")

    metrics = [
        ("new_nvars",    "nvars"),
        ("indep_sz",     "indepsz"),
        ("opt_indep_sz", "opt_isz"),
        ("irred_bin",    "irr_bin"),
        ("irred_tri",    "irr_tri"),
        ("irred_long",   "irr_long"),
        ("irred_cls",    "irr_cls"),
    ]

    union_parts = []
    for dir, ver in table_todo:
        parts = [f"replace('{dir}','out-ganak-mc','') as dirname"]
        for col, alias in metrics:
            parts.append(f"{_median_subquery(col, dir, ver, fname_like)} as med_{alias}")
            parts.append(f"{_avg_subquery(col, dir, ver, fname_like)} as avg_{alias}")
        count_sq = (f"(SELECT COUNT(*) FROM data WHERE dirname='{dir}'"
                    f" AND ganak_ver='{ver}'{fname_like})")
        parts.append(f"{count_sq} as n_inst")
        union_parts.append("SELECT " + ", ".join(parts))

    query = "\nUNION ALL\n".join(union_parts)
    if verbose:
        print(f"  Instance stats query ({len(table_todo)} rows)")
    with open("gen_table.sqlite", "w") as f:
        f.write(".mode table\n")
        f.write(query + "\n")
    os.system("sqlite3 data.sqlite3 < gen_table.sqlite")
    os.unlink("gen_table.sqlite")


def print_preproc_diffs(table_todo, fname_like, verbose=False):
    if len(table_todo) < 2:
        return

    dirs = ",".join("'" + d + "'" for d, _ in table_todo)
    vers = ",".join("'" + v + "'" for _, v in table_todo)

    con = sqlite3.connect("data.sqlite3")
    cur = con.cursor()
    cur.execute(
        f"SELECT fname FROM data"
        f" WHERE dirname IN ({dirs}) AND ganak_ver IN ({vers}){fname_like}"
        f"  AND new_nvars IS NOT NULL"
        f" GROUP BY fname"
        f" HAVING COUNT(DISTINCT new_nvars) > 1"
        f"     OR COUNT(DISTINCT indep_sz) > 1"
        f"     OR COUNT(DISTINCT opt_indep_sz) > 1"
        f"     OR COUNT(DISTINCT irred_cls) > 1"
        f" LIMIT 10"
    )
    diff_fnames = [row[0] for row in cur.fetchall()]

    if not diff_fnames:
        con.close()
        return

    title = f"Preprocessor output differences across dirs ({len(diff_fnames)} file(s) shown, up to 10)"
    print(f"\n{BLUE}{title}{RESET}")

    fnames_sql = ",".join("'" + f + "'" for f in diff_fnames)
    cur.execute(
        f"SELECT fname, replace(dirname,'out-ganak-mc',''), new_nvars, indep_sz, opt_indep_sz, irred_cls"
        f" FROM data"
        f" WHERE dirname IN ({dirs}) AND ganak_ver IN ({vers})"
        f"   AND fname IN ({fnames_sql}) AND new_nvars IS NOT NULL{fname_like}"
        f" ORDER BY fname, dirname"
    )
    rows = cur.fetchall()
    con.close()

    headers = ["fname", "dirname", "nvars", "indepsz", "opt_isz", "irred_cls"]
    str_rows = [(f, d, str(nv), str(isz), str(oisz), str(ic))
                for f, d, nv, isz, oisz, ic in rows]
    widths = [max(len(h), max((len(r[i]) for r in str_rows), default=0))
              for i, h in enumerate(headers)]
    sep = "+-" + "-+-".join("-" * w for w in widths) + "-+"
    fmt = "| " + " | ".join(f"{{:<{w}}}" for w in widths) + " |"
    print(sep)
    print(fmt.format(*headers))
    print(sep)
    prev_fname = None
    for row in str_rows:
        if prev_fname and prev_fname != row[0]:
            print(sep)
        print(fmt.format(*row))
        prev_fname = row[0]
    print(sep)


def print_two_dir_diffs(dir1, dir2, fname_like, verbose=False):
    con = sqlite3.connect("data.sqlite3")
    cur = con.cursor()
    cur.execute(
        f"SELECT a.fname,"
        f"  a.new_nvars, b.new_nvars,"
        f"  a.indep_sz, b.indep_sz,"
        f"  a.opt_indep_sz, b.opt_indep_sz,"
        f"  a.irred_cls, b.irred_cls,"
        f"  a.arjun_time, b.arjun_time"
        f" FROM data a JOIN data b ON a.fname=b.fname"
        f" WHERE a.dirname='{dir1}' AND b.dirname='{dir2}'"
        f"   AND a.new_nvars IS NOT NULL AND b.new_nvars IS NOT NULL{fname_like}"
    )
    rows = cur.fetchall()
    con.close()

    if not rows:
        return

    metric_pairs = [(1, 2), (3, 4), (7, 8)]
    if not any(r[ia] != r[ib] for r in rows for ia, ib in metric_pairs
               if r[ia] is not None and r[ib] is not None):
        return

    def diff(r, ia, ib):
        a, b = r[ia], r[ib]
        return abs((a or 0) - (b or 0))

    def top3(sort_ia, sort_ib):
        return sorted(rows, key=lambda r: diff(r, sort_ia, sort_ib), reverse=True)[:3]

    # (metric label, index_a, index_b)
    metrics = [("nvars",     1, 2),
               ("indep_sz",  3, 4),
               ("irred_cls", 7, 8)]

    # Collect top-3 per metric, preserving which metric selected each row
    sections = []
    seen = set()
    for label, ia, ib in metrics:
        section_rows = []
        for r in top3(ia, ib):
            if r[0] not in seen:
                seen.add(r[0])
            section_rows.append(r)
        sections.append((label, section_rows))

    d1s = dir1.replace("out-ganak-mc", "")
    d2s = dir2.replace("out-ganak-mc", "")

    title = f"Biggest diffs: {d1s}  vs  {d2s}"
    print(f"\n{BLUE}{title}{RESET}")
    print(f"  d1 = {d1s}")
    print(f"  d2 = {d2s}")

    def fmt_val(v):
        return "N/A" if v is None else str(v)

    def fmt_time(v):
        return "N/A" if v is None else f"{v:.1f}"

    headers = ["fname",
               "nvars-d1",  "nvars-d2",
               "indep-d1",  "indep-d2",
               "opt_i-d1",  "opt_i-d2",
               "irred-d1",  "irred-d2",
               "arjT-d1",   "arjT-d2"]

    all_str_rows = []
    for _, sec_rows in sections:
        for r in sec_rows:
            all_str_rows.append((
                r[0],
                fmt_val(r[1]),  fmt_val(r[2]),
                fmt_val(r[3]),  fmt_val(r[4]),
                fmt_val(r[5]),  fmt_val(r[6]),
                fmt_val(r[7]),  fmt_val(r[8]),
                fmt_time(r[9]), fmt_time(r[10]),
            ))

    widths = [max(len(h), max((len(row[i]) for row in all_str_rows), default=0))
              for i, h in enumerate(headers)]
    sep = "+-" + "-+-".join("-" * w for w in widths) + "-+"
    fmt_row = "| " + " | ".join(f"{{:<{w}}}" for w in widths) + " |"
    print(sep)
    print(fmt_row.format(*headers))
    print(sep)
    printed = set()
    for section_label, sec_rows in sections:
        for r in sec_rows:
            if r[0] in printed:
                continue
            printed.add(r[0])
            row_str = (
                r[0],
                fmt_val(r[1]),  fmt_val(r[2]),
                fmt_val(r[3]),  fmt_val(r[4]),
                fmt_val(r[5]),  fmt_val(r[6]),
                fmt_val(r[7]),  fmt_val(r[8]),
                fmt_time(r[9]), fmt_time(r[10]),
            )
            print(fmt_row.format(*row_str))
    print(sep)


def print_distribution(table_todo, fname_like, col, label, xscale="linear", xmin=None, xlabel=None):
    for dir, ver in table_todo:
        con = sqlite3.connect("data.sqlite3")
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

        if xmin is not None:
            values = [v for v in values if v >= xmin]
        if xscale == "log":
            values = [math.log10(v) for v in values if v > 0]

        if not values:
            continue

        title = f"{label}: {dir} [ganak_time-arjun_time >= 100s, n={len(values)}]"
        safe_dir = re.sub(r'[^a-zA-Z0-9_-]', '_', dir)
        safe_col = re.sub(r'[^a-zA-Z0-9_-]', '_', col)
        pdf_file = f"hist_{safe_col}_{safe_dir}.pdf"
        png_file = f"hist_{safe_col}_{safe_dir}.png"
        print(f"\n{BLUE}{title}{RESET}")
        print(f"  PDF: {pdf_file}  PNG: {png_file}")

        dat_file = f"hist_{safe_col}_{safe_dir}.dat"
        gp_file  = f"hist_{safe_col}_{safe_dir}.gnuplot"

        with open(dat_file, "w") as f:
            for v in values:
                f.write(f"{v}\n")

        vmin, vmax = min(values), max(values)
        binwidth = (vmax - vmin) / 20.0 if vmax > vmin else 1.0
        actual_xlabel = xlabel if xlabel is not None else col

        def gp_str(s):
            return s.replace('"', '\\"')

        with open(gp_file, "w") as f:
            for term, out in [
                ('pdfcairo size 15cm,15cm', pdf_file),
                ('pngcairo size 600,600',   png_file),
            ]:
                f.write(f'set terminal {term}\n')
                f.write(f'set output "{out}"\n')
                f.write(f'set title "{gp_str(title)}"\n')
                f.write(f'set xlabel "{gp_str(actual_xlabel)}"\n')
                f.write( 'set ylabel "count"\n')
                f.write( 'set grid\n')
                f.write( 'set key off\n')
                f.write( 'set style fill solid 0.5 border -1\n')
                f.write(f'binwidth = {binwidth}\n')
                f.write( 'set boxwidth binwidth\n')
                f.write( 'bin(x,w) = w * floor(x/w + 0.5)\n')
                f.write(f'plot "{dat_file}" using (bin($1,binwidth)):(1) smooth freq with boxes lc rgb "blue"\n\n')

        os.system(f"gnuplot {gp_file}")

        if os.path.exists(png_file):
            with open(png_file, "rb") as fh:
                img_b64 = base64.b64encode(fh.read()).decode()
            print(f"\033]1337;File=inline=1;width=600px;height=600px:{img_b64}\a")


def print_sigabrt_files(table_todo, fname_like):
    dirs = ",".join("'" + dir + "'" for dir, _ in table_todo)
    vers = ",".join("'" + ver + "'" for _, ver in table_todo)
    con = sqlite3.connect("data.sqlite3")
    cur = con.cursor()
    cur.execute(
        f"SELECT COUNT(*) FROM data WHERE dirname IN ({dirs}) AND ganak_ver IN ({vers})"
        f" AND signal=6 AND (mem_out IS NULL OR mem_out=0){fname_like}"
    )
    count = cur.fetchone()[0]
    if count == 0:
        con.close()
        return
    filter_desc = f"  |  filter: {fname_like.strip()}" if fname_like.strip() else ""
    title = f"WARNING: {count} instance(s) with sigABRT (signal=6)  |  excluding mem_out=1{filter_desc}"
    print(f"\n{BLUE}{title}{RESET}")
    cur.execute(
        f"SELECT dirname, fname, timeout_t FROM data WHERE dirname IN ({dirs}) AND ganak_ver IN ({vers})"
        f" AND signal=6 AND (mem_out IS NULL OR mem_out=0){fname_like} ORDER BY dirname, fname"
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


def print_errored_files(matched_dirs):
    if not matched_dirs:
        return
    dirs_sql = ",".join("'" + d + "'" for d in matched_dirs)
    con = sqlite3.connect("data.sqlite3")
    cur = con.cursor()
    cur.execute(
        f"SELECT COUNT(*) FROM data WHERE dirname IN ({dirs_sql}) AND errored=1"
    )
    count = cur.fetchone()[0]
    if count == 0:
        con.close()
        return
    title = f"ERROR: {count} instance(s) with ERROR or 'assertion fail' in output  |  {len(matched_dirs)} dirs matched"
    print(f"\n{RED}{title}{RESET}")
    cur.execute(
        f"SELECT dirname, fname, ganak_time FROM data WHERE dirname IN ({dirs_sql})"
        f" AND errored=1 ORDER BY dirname, fname"
    )
    rows = cur.fetchall()
    con.close()
    str_rows = [(d, f, f"{t:.2f}" if t is not None else "N/A") for d, f, t in rows]
    widths = [max(len(h), max(len(r[i]) for r in str_rows)) for i, h in enumerate(("dirname", "fname", "ganak_time"))]
    sep = RED + "+-" + "-+-".join("-" * w for w in widths) + "-+" + RESET
    fmt = RED + "| " + " | ".join(f"{{:<{w}}}" for w in widths) + " |" + RESET
    print(sep)
    print(fmt.format("dirname", "fname", "ganak_time"))
    print(sep)
    for row in str_rows:
        print(fmt.format(*row))
    print(sep)


def print_distributions(table_todo, fname_like):
    print_distribution(table_todo, fname_like, "cache_miss_rate",  "cache miss rate")
    print_distribution(table_todo, fname_like, "compsK",           "num components (K) [log10 x-axis]", xscale="log", xmin=1, xlabel="LOG compsK")
    print_distribution(table_todo, fname_like, "td_width",         "TD width", xmin=0)
    print_distribution(table_todo, fname_like, "ganak_mem_mb",     "memory usage (MB) [log10 x-axis]", xscale="log", xmin=1, xlabel="LOG mem_mb")


def scatter_plot_time_pairs(matched_dirs, fname_like, verbose=False):
    """For every pair of matched dirs, generate a gnuplot scatter plot of
    solve times (NULL -> 3600).  Writes a PDF and a PNG to disk and displays
    the PNG inline in the terminal (wezterm / iTerm2 protocol)."""
    TIMEOUT = 3600

    pairs = list(itertools.combinations(matched_dirs, 2))
    if not pairs:
        return

    con = sqlite3.connect("data.sqlite3")
    cur = con.cursor()

    for dir1, dir2 in pairs:
        query = (
            f"SELECT a.fname,"
            f" COALESCE(a.ganak_time, {TIMEOUT}),"
            f" COALESCE(b.ganak_time, {TIMEOUT})"
            f" FROM data a JOIN data b ON a.fname = b.fname"
            f" WHERE a.dirname = '{dir1}' AND b.dirname = '{dir2}'"
            f"{fname_like}"
        )
        cur.execute(query)
        rows = cur.fetchall()

        if not rows:
            if verbose:
                print(f"  scatter: no common fnames for {dir1} vs {dir2}")
            continue

        safe1 = re.sub(r'[^a-zA-Z0-9_-]', '_', dir1)
        safe2 = re.sub(r'[^a-zA-Z0-9_-]', '_', dir2)
        dat_file = f"scatter_{safe1}_vs_{safe2}.dat"
        pdf_file = f"scatter_{safe1}_vs_{safe2}.pdf"
        png_file = f"scatter_{safe1}_vs_{safe2}.png"
        gp_file  = f"scatter_{safe1}_vs_{safe2}.gnuplot"

        with open(dat_file, "w") as f:
            f.write(f"# col1={dir1}  col2={dir2}\n")
            for _fname, t1, t2 in rows:
                f.write(f"{t1}\t{t2}\n")

        # Escape double quotes for gnuplot strings
        def gp_str(s):
            return s.replace('"', '\\"')

        title  = f"Solve time: {gp_str(dir1)} vs {gp_str(dir2)}"
        xlabel = f"{gp_str(dir1)} time (s)"
        ylabel = f"{gp_str(dir2)} time (s)"

        with open(gp_file, "w") as f:
            for term, out in [
                (f'pdfcairo size 15cm,15cm', pdf_file),
                (f'pngcairo size 600,600',   png_file),
            ]:
                f.write(f'set terminal {term}\n')
                f.write(f'set output "{out}"\n')
                f.write(f'set title "{title}"\n')
                f.write(f'set xlabel "{xlabel}"\n')
                f.write(f'set ylabel "{ylabel}"\n')
                f.write( 'set logscale xy\n')
                f.write( 'set xrange [0.1:4000]\n')
                f.write( 'set yrange [0.1:4000]\n')
                f.write( 'set grid\n')
                f.write( 'set key off\n')
                f.write( 'set arrow 1 from 0.1,0.1 to 3600,3600 nohead lc rgb "gray50" lw 1\n')
                f.write(f'plot "{dat_file}" using 1:2 with points pt 7 ps 0.5 lc rgb "blue" notitle\n')
                f.write( 'unset arrow 1\n\n')

        console_title = f"Scatter solve time: {dir1} vs {dir2} (n={len(rows)})"
        print(f"\n{BLUE}{console_title}{RESET}")
        print(f"  PDF: {pdf_file}  PNG: {png_file}")

        os.system(f"gnuplot {gp_file}")

        # Display PNG inline (wezterm / iTerm2 inline-image protocol)
        if os.path.exists(png_file):
            with open(png_file, "rb") as fh:
                img_b64 = base64.b64encode(fh.read()).decode()
            print(f"\033]1337;File=inline=1;width=600px;height=600px:{img_b64}\a")

    con.close()


def generate_gnuplot(fname2_s, verbose=False):
    gnuplotfn = "cdf.gnuplot"
    pdf_file = "cdf.pdf"
    png_file = "cdf.png"
    if verbose:
        print(f"Writing gnuplot script to {gnuplotfn} with {len(fname2_s)} data series")

    def gp_str(s):
        return s.replace('"', '\\"')

    def plot_lines():
        towrite = ""
        for fn, call, ver, _, dir in fname2_s:
            if True:
                call = gnuplot_name_cleanup(call)
                dir  = gnuplot_name_cleanup(dir)
                ver  = gnuplot_name_cleanup(ver)
                oneline = "\""+fn+"\" u 2:1 with linespoints  title \""+ver+"-"+dir+"-"+call+"\""
                towrite += oneline
                towrite += ",\\\n"
        return towrite[:(len(towrite)-4)]

    with open(gnuplotfn, "w") as f:
        for term, out in [
            ('pdfcairo size 15cm,15cm', pdf_file),
            ('pngcairo size 600,600',   png_file),
        ]:
            f.write(f'set terminal {term}\n')
            f.write(f'set output "{out}"\n')
            f.write('set title "Counter ganak"\n')
            f.write('set key bottom right\n')
            f.write('set logscale x\n')
            f.write('unset logscale y\n')
            f.write('set ylabel "Instances counted"\n')
            f.write('set xlabel "Time (s)"\n')
            f.write('set grid\n')
            f.write('plot [0.1:3600][0.1:]\\\n')
            f.write(plot_lines())
            f.write('\n\n')
    return gnuplotfn, pdf_file, png_file


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
conn = sqlite3.connect('data.sqlite3')
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
conn = sqlite3.connect('data.sqlite3')

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
    "out-ganak-mccomp2324-1193808-0", # fixing counting bug, adding new cube extend system, fixing cube counting (and maybe effectiveness)
    # "out-ganak-mccomp2324-1229753-0", # lots of bug fixes, beauty changes with Claude, etc
    # "out-ganak-mccomp2324-1231407-0", # the same as above but without (most) of the Claude improvements
    # "out-ganak-mccomp2324-1247484-0", ## new shrinking, fixing arjun SLOW_DEBUG, improved propagation idea from CaDiCaL, improved 3-tier clause database. Also undoing only 2 particular Claude changes (propagation and --prob 0 non-zeroing of data)
    # "out-ganak-mccomp2324-1250247-", # CMS cleanup, oracle improvements, fix parsing issue (?) of CNF header -- weird + one of them is a binary with: reason-side bumping and lbd update, evsids
    "out-ganak-mccomp2324-1256426-0", # fixing oracle, mostly, and also header parsing more lax
    "out-ganak-mccomp2324-1261017-0", # Different order in Arjun
]
only_dirs = [
     "mei-march-2026-1239767-1", # gpmc
     # "mei-march-2026-1239767-0", # ganak old
     # "mei-march-2026-1269673-0", # ganak release, but WRONG SED
     "mei-march-2026-1274973-0", # ganak release, new SED
]

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
    parser = argparse.ArgumentParser(description="Generate CDF plots and tables for ganak benchmark data")
    parser.add_argument("--verbose", "-v", action="store_true", help="Print progress information")
    parser.add_argument("--full", action="store_true", help="Print full summary table (default: compact)")

    parser.add_argument("--fname", nargs="+", metavar="PATTERN", default=[],
                        help="Filter by fname pattern(s), e.g. --fname '%%track1%%' '%%track3%%'")
    args = parser.parse_args()

    if args.fname:
        clauses = " or ".join(f"fname like '{p}'" for p in args.fname)
        fname_like = f" and ({clauses}) "
    else:
        fname_like = ""

    matched_dirs = get_matching_dirs(only_dirs)
    if args.verbose:
        print(f"Found {len(versions)} versions in database")
        print(f"Matched {len(matched_dirs)} dirs from only_dirs prefixes")
        print("Building CSV data...")
    scatter_plot_time_pairs(matched_dirs, fname_like, args.verbose)
    fname2_s, table_todo = build_csv_data(todo, matched_dirs, only_calls, not_calls, not_versions, fname_like, args.verbose)

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
    print_errored_files(matched_dirs)

    if args.verbose:
        print("Printing median tables...")
    print_median_tables(table_todo, fname_like, args.verbose)
    print_instance_stats_table(table_todo, fname_like, args.verbose)
    print_preproc_diffs(table_todo, fname_like, args.verbose)
    unique_dirs = list(dict.fromkeys(d for d, _ in table_todo))
    for dir1, dir2 in itertools.combinations(unique_dirs, 2):
        print_two_dir_diffs(dir1, dir2, fname_like, args.verbose)

    if args.verbose:
        print("Generating gnuplot script...")
    gnuplotfn, pdf_file, png_file = generate_gnuplot(fname2_s, args.verbose)

    dirs = ",".join("'" + dir + "'" for dir, _ in table_todo)
    create_notebook(dirs)

    for path in [pdf_file, png_file]:
        if os.path.exists(path):
            os.unlink(path)
    os.system(f"gnuplot {gnuplotfn}")

    console_title = f"CDF: instances counted vs. solve time"
    print(f"\n{BLUE}{console_title}{RESET}")
    print(f"  PDF: {pdf_file}  PNG: {png_file}")

    if os.path.exists(png_file):
        with open(png_file, "rb") as fh:
            img_b64 = base64.b64encode(fh.read()).decode()
        print(f"\033]1337;File=inline=1;width=600px;height=600px:{img_b64}\a")


if __name__ == "__main__":
    main()
