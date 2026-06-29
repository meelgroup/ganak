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
GREEN  = "\033[92m"
RED    = "\033[91m"
RESET = "\033[0m"

TMP_DIR = "tmp"


def convert_to_cdf(fname, fname2):
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
            fname = f"{TMP_DIR}/run-"+dir+".csv"
            gencsv = f"{TMP_DIR}/gencsv.sqlite"
            with open(gencsv, "w") as f:
                f.write(".headers off\n")
                f.write(".mode csv\n")
                f.write(".output "+fname+"\n")
                f.write("select ganak_time from data where dirname='"+dir+"' and ganak_ver='"+ver+"'\n and ganak_time is not NULL "+fname_like)
            os.system(f"sqlite3 data.sqlite3 < {gencsv}")
            os.unlink(gencsv)

            fname2 = fname + ".gnuplotdata"
            num_solved = convert_to_cdf(fname, fname2)
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
        gen_table = f"{TMP_DIR}/gen_table.sqlite"
        with open(gen_table, "w") as f:
            f.write(".mode table\n")
            # f.write(".mode colum\n")
            # f.write(".headers off\n")
            query = (f"select\n        {select_clause}\n"
                     f"        from data where dirname IN ({dirs}) and ganak_ver IN ({vers})"
                     f" {fname_like} {counted_req}group by dirname order by solved asc")
            if verbose:
                print(f"  Summary query: {query[:120]}...")
            f.write(query)
        os.system(f"sqlite3 data.sqlite3 < {gen_table}")
        os.unlink(gen_table)


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
    gen_table = f"{TMP_DIR}/gen_table.sqlite"
    with open(gen_table, "w") as f:
        f.write(".mode table\n")
        f.write(query + "\n")
    os.system(f"sqlite3 data.sqlite3 < {gen_table}")
    os.unlink(gen_table)


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

    for min_nvars, label in [(None, "all"), (20, "nvars>=20")]:
        nvars_filter = f" and new_nvars >= {min_nvars}" if min_nvars is not None else ""
        fl = fname_like + nvars_filter

        union_parts = []
        for dir, ver in table_todo:
            parts = [f"replace('{dir}','out-ganak-mc','') as dirname"]
            for col, alias in metrics:
                parts.append(f"{_median_subquery(col, dir, ver, fl)} as med_{alias}")
                parts.append(f"{_avg_subquery(col, dir, ver, fl)} as avg_{alias}")
            count_sq = (f"(SELECT COUNT(*) FROM data WHERE dirname='{dir}'"
                        f" AND ganak_ver='{ver}'{fl})")
            parts.append(f"{count_sq} as n_inst")
            union_parts.append("SELECT " + ", ".join(parts))

        query = "\nUNION ALL\n".join(union_parts)
        print(f"\n  [{label}]")
        if verbose:
            print(f"  Instance stats query ({len(table_todo)} rows)")
        gen_table = f"{TMP_DIR}/gen_table.sqlite"
        with open(gen_table, "w") as f:
            f.write(".mode table\n")
            f.write(query + "\n")
        os.system(f"sqlite3 data.sqlite3 < {gen_table}")
        os.unlink(gen_table)


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
        f"SELECT fname, replace(dirname,'out-ganak-mc',''), new_nvars, indep_sz, opt_indep_sz, irred_cls, arjun_time, ganak_time"
        f" FROM data"
        f" WHERE dirname IN ({dirs}) AND ganak_ver IN ({vers})"
        f"   AND fname IN ({fnames_sql}) AND new_nvars IS NOT NULL{fname_like}"
        f" ORDER BY fname, dirname"
    )
    rows = cur.fetchall()
    con.close()

    headers = ["fname", "dirname", "nvars", "indepsz", "opt_isz", "irred_cls", "arjun_t", "ganak_t"]
    str_rows = [(f, d, str(nv), str(isz), str(oisz), str(ic),
                 f"{at:.2f}" if at is not None else "NULL",
                 f"{gt:.2f}" if gt is not None else "NULL")
                for f, d, nv, isz, oisz, ic, at, gt in rows]
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


def print_solved_only_diffs(dir1, dir2, fname_like, verbose=False):
    """For the pair (dir1, dir2), print two tables:
      - files solved by dir1 but NOT by dir2
      - files solved by dir2 but NOT by dir1
    'Solved' means ganak_time IS NOT NULL.
    """
    con = sqlite3.connect("data.sqlite3")
    cur = con.cursor()
    cur.execute(
        f"SELECT fname, dirname, new_nvars, indep_sz, opt_indep_sz,"
        f"       irred_cls, ganak_time, arjun_time, td_width"
        f" FROM data WHERE dirname IN ('{dir1}','{dir2}'){fname_like}"
    )
    per_dir = {dir1: {}, dir2: {}}
    for fname, dirname, nvars, isz, oisz, irred, gt, at, td in cur.fetchall():
        per_dir[dirname][fname] = (nvars, isz, oisz, irred, gt, at, td)
    con.close()

    def solved(d, f):
        row = per_dir[d].get(f)
        return row is not None and row[4] is not None

    all_fnames = set(per_dir[dir1]) | set(per_dir[dir2])
    only_in_1 = sorted(f for f in all_fnames if solved(dir1, f) and not solved(dir2, f))
    only_in_2 = sorted(f for f in all_fnames if solved(dir2, f) and not solved(dir1, f))

    d1s = dir1.replace("out-ganak-mc", "")
    d2s = dir2.replace("out-ganak-mc", "")

    def fmt_val(v):
        return "N/A" if v is None else str(v)
    def fmt_time(v):
        return "N/A" if v is None else f"{v:.1f}"

    headers = ["fname",
               "nvars-d1", "nvars-d2",
               "indep-d1", "indep-d2",
               "opt_i-d1", "opt_i-d2",
               "irred-d1", "irred-d2",
               "gnkT-d1",  "gnkT-d2",
               "arjT-d1",  "arjT-d2",
               "td-d1",    "td-d2"]

    def build_row(f):
        a = per_dir[dir1].get(f, (None,)*7)
        b = per_dir[dir2].get(f, (None,)*7)
        return (f,
                fmt_val(a[0]),  fmt_val(b[0]),
                fmt_val(a[1]),  fmt_val(b[1]),
                fmt_val(a[2]),  fmt_val(b[2]),
                fmt_val(a[3]),  fmt_val(b[3]),
                fmt_time(a[4]), fmt_time(b[4]),
                fmt_time(a[5]), fmt_time(b[5]),
                fmt_val(a[6]),  fmt_val(b[6]))

    def print_section(label, fnames):
        str_rows = [build_row(f) for f in fnames]
        widths = [max(len(h), max((len(row[i]) for row in str_rows), default=0))
                  for i, h in enumerate(headers)]
        sep = "+-" + "-+-".join("-" * w for w in widths) + "-+"
        fmt_row = "| " + " | ".join(f"{{:<{w}}}" for w in widths) + " |"
        print(f"\n{BLUE}{label}{RESET}")
        print(f"  d1 = {d1s}")
        print(f"  d2 = {d2s}")
        print(sep)
        print(fmt_row.format(*headers))
        print(sep)
        if not str_rows:
            print(f"{GREEN}EMPTY{RESET}")
        else:
            for row in str_rows:
                print(fmt_row.format(*row))
        print(sep)

    print_section(
        f"Solved only by d1 ({d1s}), not by d2 ({d2s})  [{len(only_in_1)} files]",
        only_in_1)
    print_section(
        f"Solved only by d2 ({d2s}), not by d1 ({d1s})  [{len(only_in_2)} files]",
        only_in_2)


def print_solution_mismatches(dir1, dir2, fname_like, verbose=False):
    """For the pair (dir1, dir2), print files where BOTH solved but the
    log10-estimate of the model count differs. Table is all red."""
    con = sqlite3.connect("data.sqlite3")
    cur = con.cursor()
    cur.execute(
        f"SELECT a.fname,"
        f"  a.new_nvars, b.new_nvars,"
        f"  a.indep_sz, b.indep_sz,"
        f"  a.opt_indep_sz, b.opt_indep_sz,"
        f"  a.irred_cls, b.irred_cls,"
        f"  a.ganak_time, b.ganak_time,"
        f"  a.arjun_time, b.arjun_time,"
        f"  a.mc_log10, b.mc_log10"
        f" FROM data a JOIN data b ON a.fname=b.fname"
        f" WHERE a.dirname='{dir1}' AND b.dirname='{dir2}'"
        f"   AND a.ganak_time IS NOT NULL AND b.ganak_time IS NOT NULL"
        f"   AND a.mc_log10 IS NOT NULL AND b.mc_log10 IS NOT NULL{fname_like}"
    )
    rows = cur.fetchall()
    con.close()

    # Compare with a small tolerance to avoid floating-point noise.
    def mismatched(a, b):
        return abs(a - b) > 1e-6 * max(1.0, abs(a), abs(b))

    bad = [r for r in rows if mismatched(r[13], r[14])]
    d1s = dir1.replace("out-ganak-mc", "")
    d2s = dir2.replace("out-ganak-mc", "")

    def fmt_val(v):
        return "N/A" if v is None else str(v)
    def fmt_time(v):
        return "N/A" if v is None else f"{v:.1f}"
    def fmt_log(v):
        return "N/A" if v is None else f"{v:.6g}"

    headers = ["fname",
               "nvars-d1", "nvars-d2",
               "indep-d1", "indep-d2",
               "opt_i-d1", "opt_i-d2",
               "irred-d1", "irred-d2",
               "gnkT-d1",  "gnkT-d2",
               "arjT-d1",  "arjT-d2",
               "log10-d1", "log10-d2"]

    str_rows = []
    for r in bad:
        str_rows.append((
            r[0],
            fmt_val(r[1]),  fmt_val(r[2]),
            fmt_val(r[3]),  fmt_val(r[4]),
            fmt_val(r[5]),  fmt_val(r[6]),
            fmt_val(r[7]),  fmt_val(r[8]),
            fmt_time(r[9]), fmt_time(r[10]),
            fmt_time(r[11]), fmt_time(r[12]),
            fmt_log(r[13]), fmt_log(r[14]),
        ))

    widths = [max(len(h), max((len(row[i]) for row in str_rows), default=0))
              for i, h in enumerate(headers)]
    sep = "+-" + "-+-".join("-" * w for w in widths) + "-+"
    fmt_row = "| " + " | ".join(f"{{:<{w}}}" for w in widths) + " |"

    title = f"SOLUTION MISMATCH (both solved, different log10): {d1s} vs {d2s}  [{len(bad)} files]"
    print(f"\n{RED}{title}{RESET}")
    print(f"{RED}  d1 = {d1s}{RESET}")
    print(f"{RED}  d2 = {d2s}{RESET}")
    print(f"{RED}{sep}{RESET}")
    print(f"{RED}{fmt_row.format(*headers)}{RESET}")
    print(f"{RED}{sep}{RESET}")
    if not str_rows:
        print(f"{GREEN}EMPTY{RESET}")
    else:
        for row in str_rows:
            print(f"{RED}{fmt_row.format(*row)}{RESET}")
    print(f"{RED}{sep}{RESET}")


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
        pdf_file = f"{TMP_DIR}/hist_{safe_col}_{safe_dir}.pdf"
        png_file = f"{TMP_DIR}/hist_{safe_col}_{safe_dir}.png"
        print(f"\n{BLUE}{title}{RESET}")
        print(f"  PDF: {pdf_file}  PNG: {png_file}")

        dat_file = f"{TMP_DIR}/hist_{safe_col}_{safe_dir}.dat"
        gp_file  = f"{TMP_DIR}/hist_{safe_col}_{safe_dir}.gnuplot"

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
                ('pdfcairo size 52cm,40cm background "#d0d0d0"', pdf_file),
                ('pngcairo size 600,600 background "#d0d0d0"',   png_file),
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
        _display_png(png_file)


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


def _preproc_has_data(matched_dirs):
    """Return True if the preproc table exists and has data for any of the matched dirs."""
    con = sqlite3.connect("data.sqlite3")
    cur = con.cursor()
    try:
        cur.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='preproc'")
        if not cur.fetchone():
            return False
        if not matched_dirs:
            return False
        dirs_sql = ",".join("'" + d + "'" for d in matched_dirs)
        cur.execute(f"SELECT COUNT(*) FROM preproc WHERE dirname IN ({dirs_sql})")
        count = cur.fetchone()[0]
        return count > 0
    except Exception:
        return False
    finally:
        con.close()


def _preproc_step_stats(con, dirs_sql, where_extra="", per_cnf=False):
    """Fetch all preproc delta rows and compute per-step stats including medians in Python.

    If per_cnf=True, first sum all occurrences within each (dirname, fname) before computing
    the median across CNFs — giving the total effect per CNF per step.

    Returns dict: name -> {n, n_active, med_lits, med_bins, med_cls, med_vars,
                            med_step_num}
    """
    from collections import defaultdict
    cur = con.cursor()
    cur.execute(
        f"SELECT name, dirname, fname, delta_irred_bins, delta_irred_long_cls,"
        f" delta_irred_long_lits, delta_free_vars, step_num"
        f" FROM preproc WHERE dirname IN ({dirs_sql}) AND depth = 0{where_extra}"
    )

    if per_cnf:
        # Accumulate per (name, dirname, fname), then collect per name
        # cnf_totals[name][(dirname, fname)] = {'bins': int, 'cls': int, 'lits': int, 'vars': int, 'snums': list}
        cnf_totals: dict = {}
        for name, dirname, fname, d_bins, d_cls, d_lits, d_vars, snum in cur.fetchall():
            key = (dirname, fname)
            if name not in cnf_totals:
                cnf_totals[name] = {}
            if key not in cnf_totals[name]:
                cnf_totals[name][key] = {'bins': 0, 'cls': 0, 'lits': 0, 'vars': 0, 'snums': []}
            t = cnf_totals[name][key]
            t['bins'] += d_bins if d_bins is not None else 0
            t['cls']  += d_cls  if d_cls  is not None else 0
            t['lits'] += d_lits if d_lits is not None else 0
            t['vars'] += d_vars if d_vars is not None else 0
            if snum is not None:
                t['snums'].append(snum)
        groups = defaultdict(lambda: {'bins': [], 'cls': [], 'lits': [], 'vars': [], 'snums': []})
        for name, cnf_dict in cnf_totals.items():
            g = groups[name]
            for t in cnf_dict.values():
                g['bins'].append(t['bins'])
                g['cls'].append(t['cls'])
                g['lits'].append(t['lits'])
                g['vars'].append(t['vars'])
                g['snums'].extend(t['snums'])
    else:
        groups = defaultdict(lambda: {'bins': [], 'cls': [], 'lits': [], 'vars': [], 'snums': []})
        for name, _, _, d_bins, d_cls, d_lits, d_vars, snum in cur.fetchall():
            g = groups[name]
            g['bins'].append(d_bins if d_bins is not None else 0)
            g['cls'].append(d_cls  if d_cls  is not None else 0)
            g['lits'].append(d_lits if d_lits is not None else 0)
            g['vars'].append(d_vars if d_vars is not None else 0)
            if snum is not None:
                g['snums'].append(snum)

    def median(vals):
        s = sorted(vals)
        return s[len(s) // 2] if s else 0

    result = {}
    for name, g in groups.items():
        lits = g['lits']
        vars_ = g['vars']
        result[name] = {
            'n':            len(lits),
            'med_bins':     median(g['bins']),
            'med_cls':      median(g['cls']),
            'med_lits':     median(lits),
            'med_vars':     median(vars_),
            'med_step_num': median(g['snums']) if g['snums'] else 0,
        }
    return result


def _print_table(headers, str_rows):
    widths = [max(len(h), max((len(r[i]) for r in str_rows), default=0))
              for i, h in enumerate(headers)]
    sep = "+-" + "-+-".join("-" * w for w in widths) + "-+"
    fmt = "| " + " | ".join(f"{{:<{w}}}" for w in widths) + " |"
    print(sep)
    print(fmt.format(*headers))
    print(sep)
    for row in str_rows:
        print(fmt.format(*row))
    print(sep)


def _preproc_file_label(matched_dirs):
    """Return a filesystem-safe label for chart filenames, derived from the dirname."""
    if len(matched_dirs) == 1:
        return re.sub(r'[^a-zA-Z0-9_-]', '_', matched_dirs[0])[-50:]
    return "combined"


def _preproc_has_step_time(matched_dirs):
    """Return True if step_time column exists and has data for matched dirs."""
    try:
        con = sqlite3.connect("data.sqlite3")
        cur = con.cursor()
        dirs_sql = ",".join("'" + d + "'" for d in matched_dirs)
        cur.execute(f"SELECT COUNT(*) FROM preproc WHERE dirname IN ({dirs_sql}) AND step_time IS NOT NULL")
        count = cur.fetchone()[0]
        con.close()
        return count > 0
    except Exception:
        return False


def print_preproc_delta_table(matched_dirs, verbose=False):
    """Per-step total contribution: SUM of each delta across all CNFs, sorted by lits impact."""
    if not _preproc_has_data(matched_dirs):
        return

    has_time = _preproc_has_step_time(matched_dirs)
    dirs_sql = ",".join("'" + d + "'" for d in matched_dirs)
    con = sqlite3.connect("data.sqlite3")
    cur = con.cursor()
    cur.execute(f"""
        SELECT name,
               COUNT(*)                                                       as n_calls,
               SUM(delta_irred_long_lits)                                     as sum_d_lits,
               SUM(delta_irred_long_cls)                                      as sum_d_cls,
               SUM(delta_irred_bins)                                          as sum_d_bins,
               SUM(delta_free_vars)                                           as sum_d_vars,
               ROUND(100.0 * SUM(CASE WHEN delta_irred_long_lits < 0
                                 THEN 1 ELSE 0 END) / COUNT(*), 0)            as pct_invoc_red_lits,
               ROUND(100.0 * SUM(CASE WHEN delta_free_vars < 0
                                 THEN 1 ELSE 0 END) / COUNT(*), 0)            as pct_invoc_red_vars,
               ROUND(SUM(step_time), 2)                                       as total_step_s
        FROM preproc WHERE dirname IN ({dirs_sql}) AND depth = 0
        GROUP BY name
        ORDER BY sum_d_lits ASC
    """)
    rows = cur.fetchall()
    con.close()

    if not rows:
        return

    title = "Preprocessing step total contribution (sum across all invocations, sorted by lits removed)"
    print(f"\n{BLUE}{title}{RESET}")
    print("  (n_calls = total invocations; d_lits/d_cls/d_bins/d_fvars = raw delta,"
          " negative = removed; %calls_lit = % of calls that reduced lits;"
          " %calls_vars = % of calls that removed a var; total_s = wall time)")

    if has_time:
        headers = ["step", "n_calls", "d_lits", "d_long_cls", "d_bin_cls", "d_fvars", "%calls_lit", "%calls_vars", "total_s"]
        str_rows = [
            (str(r[0]), str(r[1]), str(int(r[2] or 0)), str(int(r[3] or 0)),
             str(int(r[4] or 0)), str(int(r[5] or 0)), f"{int(r[6] or 0)}%",
             f"{int(r[7] or 0)}%",
             f"{r[8]:.1f}" if r[8] is not None else "N/A")
            for r in rows
        ]
    else:
        headers = ["step", "n_calls", "d_lits", "d_long_cls", "d_bin_cls", "d_fvars", "%calls_lit", "%calls_vars"]
        str_rows = [
            (str(r[0]), str(r[1]), str(int(r[2] or 0)), str(int(r[3] or 0)),
             str(int(r[4] or 0)), str(int(r[5] or 0)), f"{int(r[6] or 0)}%",
             f"{int(r[7] or 0)}%")
            for r in rows
        ]
    _print_table(headers, str_rows)


def print_preproc_step_efficiency(matched_dirs):
    """Step efficiency: total lits/vars removed (SUM), % of invocations that do anything,
    and average removal per active invocation. Sorted by total lits removed."""
    if not _preproc_has_data(matched_dirs):
        return

    has_time = _preproc_has_step_time(matched_dirs)
    dirs_sql = ",".join("'" + d + "'" for d in matched_dirs)
    con = sqlite3.connect("data.sqlite3")
    cur = con.cursor()
    cur.execute(f"""
        SELECT name,
               COUNT(*)                                                              AS n_calls,
               -SUM(delta_irred_long_lits)                                           AS lits_rmvd,
               -SUM(delta_irred_long_cls)                                            AS cls_rmvd,
               -SUM(delta_free_vars)                                                 AS fvars_rmvd,
               SUM(CASE WHEN delta_irred_long_lits < 0 THEN 1 ELSE 0 END)           AS calls_lit,
               SUM(CASE WHEN delta_irred_long_cls < 0 THEN 1 ELSE 0 END)            AS calls_cls,
               SUM(CASE WHEN delta_free_vars < 0 THEN 1 ELSE 0 END)                 AS calls_var,
               ROUND(SUM(step_time), 2)                                              AS total_s
        FROM preproc
        WHERE dirname IN ({dirs_sql}) AND depth = 0
        GROUP BY name
        HAVING lits_rmvd > 0 OR fvars_rmvd > 0
        ORDER BY lits_rmvd DESC
    """)
    rows = cur.fetchall()
    con.close()

    if not rows:
        return

    title = "Preprocessing step efficiency (sum across all invocations)"
    print(f"\n{BLUE}{title}{RESET}")
    print("  (lits_rmvd = irred long lits removed; fvars_rmvd = free vars freed;"
          " %lits_act/%cls_act/%vars_act = % of calls with nonzero removal;"
          " lits/act, cls/act, vars/act = avg removed per active call)")

    if has_time:
        headers = ["step", "n_calls", "lits_rmvd", "fvars_rmvd", "%lits_act", "%cls_act", "%vars_act",
                   "lits/act", "cls/act", "vars/act", "total_s"]
        str_rows = []
        for name, n_calls, lits_rmvd, cls_rmvd, fvars_rmvd, calls_lit, calls_cls, calls_var, total_s in rows:
            avg_lits = int(lits_rmvd / calls_lit) if calls_lit else 0
            avg_cls  = int(cls_rmvd  / calls_cls)  if calls_cls  else 0
            avg_vars = int(fvars_rmvd / calls_var) if calls_var else 0
            pct_lits = f"{100*calls_lit//n_calls}%" if n_calls else "0%"
            pct_cls  = f"{100*calls_cls//n_calls}%" if n_calls else "0%"
            pct_vars = f"{100*calls_var//n_calls}%" if n_calls else "0%"
            str_rows.append((name, str(n_calls), f"{lits_rmvd:,}", f"{fvars_rmvd:,}",
                             pct_lits, pct_cls, pct_vars,
                             f"{avg_lits:,}", f"{avg_cls:,}", f"{avg_vars:,}",
                             f"{total_s:.1f}" if total_s is not None else "N/A"))
    else:
        headers = ["step", "n_calls", "lits_rmvd", "fvars_rmvd", "%lits_act", "%cls_act", "%vars_act",
                   "lits/act", "cls/act", "vars/act"]
        str_rows = []
        for name, n_calls, lits_rmvd, cls_rmvd, fvars_rmvd, calls_lit, calls_cls, calls_var, _ in rows:
            avg_lits = int(lits_rmvd / calls_lit) if calls_lit else 0
            avg_cls  = int(cls_rmvd  / calls_cls)  if calls_cls  else 0
            avg_vars = int(fvars_rmvd / calls_var) if calls_var else 0
            pct_lits = f"{100*calls_lit//n_calls}%" if n_calls else "0%"
            pct_cls  = f"{100*calls_cls//n_calls}%" if n_calls else "0%"
            pct_vars = f"{100*calls_var//n_calls}%" if n_calls else "0%"
            str_rows.append((name, str(n_calls), f"{lits_rmvd:,}", f"{fvars_rmvd:,}",
                             pct_lits, pct_cls, pct_vars,
                             f"{avg_lits:,}", f"{avg_cls:,}", f"{avg_vars:,}"))
    _print_table(headers, str_rows)


def print_preproc_time_breakdown(matched_dirs):
    """Time breakdown: total seconds and % of preprocessing time per step, sorted by time."""
    if not _preproc_has_data(matched_dirs) or not _preproc_has_step_time(matched_dirs):
        return

    dirs_sql = ",".join("'" + d + "'" for d in matched_dirs)
    con = sqlite3.connect("data.sqlite3")
    cur = con.cursor()
    cur.execute(f"""
        SELECT name,
               COUNT(*)                                                               AS n_calls,
               ROUND(SUM(step_time), 2)                                               AS total_s,
               ROUND(AVG(step_time) * 1000, 1)                                        AS avg_ms,
               -SUM(delta_irred_long_lits)                                             AS lits_rmvd,
               -SUM(delta_irred_long_cls)                                              AS cls_rmvd,
               -SUM(delta_free_vars)                                                   AS fvars_rmvd,
               SUM(IFNULL(delta_elimed_vars, 0))                                       AS elim_vars,
               SUM(IFNULL(delta_units, 0))                                             AS units_found
        FROM preproc
        WHERE dirname IN ({dirs_sql}) AND depth = 0 AND step_time IS NOT NULL
        GROUP BY name
        ORDER BY total_s DESC
    """)
    rows = cur.fetchall()
    con.close()

    if not rows:
        return

    total_t = sum(r[2] for r in rows if r[2] is not None)

    title = "Preprocessing time breakdown per step (sorted by time)"
    print(f"\n{BLUE}{title}{RESET}")
    print(f"  Total preprocessing time: {total_t:.1f}s across {len(matched_dirs)} dirs")
    print("  (n_calls = total invocations; lits_rmvd = irred long lits removed;"
          " cls_rmvd = irred long clauses removed; fvars_rmvd = free vars freed;"
          " elim_vars = vars eliminated by BVE;"
          " units = unit clauses found; lits/s, cls/s = removed per second of step time)")

    headers = ["step", "n_calls", "total_s", "%time", "avg_ms", "lits_rmvd", "cls_rmvd",
               "fvars_rmvd", "elim_vars", "units", "lits/s", "cls/s"]
    str_rows = []
    for name, n_calls, total_s, avg_ms, lits_rmvd, cls_rmvd, fvars_rmvd, elim_vars, units_found in rows:
        pct = f"{100.0 * total_s / total_t:.1f}%" if total_t else "0%"
        lits_per_s = int(max(lits_rmvd, 0) / total_s) if total_s and total_s > 0 else 0
        cls_per_s  = int(max(cls_rmvd,  0) / total_s) if total_s and total_s > 0 else 0
        str_rows.append((
            name,
            str(n_calls),
            f"{total_s:.1f}",
            pct,
            f"{avg_ms:.1f}",
            f"{max(lits_rmvd, 0):,}",
            f"{max(cls_rmvd,  0):,}",
            f"{max(fvars_rmvd, 0):,}",
            f"{elim_vars:,}" if elim_vars else "0",
            f"{units_found:,}" if units_found else "0",
            f"{lits_per_s:,}",
            f"{cls_per_s:,}",
        ))
    _print_table(headers, str_rows)


def print_preproc_noop_waste(matched_dirs):
    """No-op waste: time spent in steps that do absolutely nothing. Potential optimization target."""
    if not _preproc_has_data(matched_dirs) or not _preproc_has_step_time(matched_dirs):
        return

    dirs_sql = ",".join("'" + d + "'" for d in matched_dirs)
    con = sqlite3.connect("data.sqlite3")
    cur = con.cursor()
    cur.execute(f"""
        SELECT name,
               COUNT(*)                                                                     AS n_total,
               SUM(CASE WHEN delta_irred_long_lits = 0 AND delta_free_vars = 0
                             AND IFNULL(delta_elimed_vars, 0) = 0
                             AND IFNULL(delta_replaced_vars, 0) = 0
                             AND IFNULL(delta_units, 0) = 0
                        THEN 1 ELSE 0 END)                                                  AS n_noop,
               ROUND(SUM(CASE WHEN delta_irred_long_lits = 0 AND delta_free_vars = 0
                             AND IFNULL(delta_elimed_vars, 0) = 0
                             AND IFNULL(delta_replaced_vars, 0) = 0
                             AND IFNULL(delta_units, 0) = 0
                        THEN step_time ELSE 0 END), 2)                                      AS wasted_s,
               ROUND(SUM(step_time), 2)                                                     AS total_s
        FROM preproc
        WHERE dirname IN ({dirs_sql}) AND depth = 0 AND step_time IS NOT NULL
        GROUP BY name
        HAVING n_noop > 0
        ORDER BY wasted_s DESC
    """)
    rows = cur.fetchall()
    con.close()

    if not rows:
        return

    total_wasted = sum(r[3] for r in rows if r[3] is not None)  # r[3] = wasted_s
    title = "No-op waste: time in steps that change nothing (potential optimization targets)"
    print(f"\n{BLUE}{title}{RESET}")
    print(f"  Total wasted time: {total_wasted:.1f}s")

    headers = ["step", "n_total", "n_noop", "%noop", "wasted_s", "total_s", "%wasted_of_step"]
    str_rows = []
    for name, n_total, n_noop, wasted_s, total_s in rows:
        pct_noop = f"{100 * n_noop // n_total}%" if n_total else "0%"
        pct_wasted = f"{100.0 * wasted_s / total_s:.0f}%" if total_s else "0%"
        str_rows.append((name, str(n_total), str(n_noop), pct_noop,
                         f"{wasted_s:.1f}", f"{total_s:.1f}", pct_wasted))
    _print_table(headers, str_rows)


def print_preproc_extended_vars(matched_dirs):
    """Extended var analysis: elim_vars, replaced_vars, units found per step."""
    if not _preproc_has_data(matched_dirs):
        return

    dirs_sql = ",".join("'" + d + "'" for d in matched_dirs)
    con = sqlite3.connect("data.sqlite3")
    cur = con.cursor()
    cur.execute(f"""
        SELECT name,
               COUNT(*)                                   AS n_active,
               SUM(IFNULL(delta_elimed_vars, 0))           AS sum_elim,
               SUM(IFNULL(delta_replaced_vars, 0))         AS sum_repl,
               SUM(IFNULL(delta_units, 0))                 AS sum_units,
               ROUND(SUM(step_time), 2)                    AS total_s
        FROM preproc
        WHERE dirname IN ({dirs_sql}) AND depth = 0
          AND (IFNULL(delta_elimed_vars,0) != 0
            OR IFNULL(delta_replaced_vars,0) != 0
            OR IFNULL(delta_units,0) != 0)
        GROUP BY name
        ORDER BY (sum_elim + sum_repl + sum_units) DESC
    """)
    rows = cur.fetchall()
    con.close()

    if not rows:
        return

    has_time = _preproc_has_step_time(matched_dirs)
    title = "Extended variable analysis: elim/replaced/unit propagations per step"
    print(f"\n{BLUE}{title}{RESET}")

    print("  (n_active = invocations where step changed at least one of elim/repl/units; "
          "total_s = time across those active invocations only)")
    if has_time:
        headers = ["step", "n_active", "elim_vars", "repl_vars", "units_found", "total_s"]
        str_rows = [(name, str(n), f"{e:,}", f"{r:,}", f"{u:,}",
                     f"{t:.1f}" if t is not None else "N/A")
                    for name, n, e, r, u, t in rows]
    else:
        headers = ["step", "n_active", "elim_vars", "repl_vars", "units_found"]
        str_rows = [(name, str(n), f"{e:,}", f"{r:,}", f"{u:,}")
                    for name, n, e, r, u, _ in rows]
    _print_table(headers, str_rows)


def print_step_predecessor_effectiveness(matched_dirs, step="must-scc-vrepl"):
    """For a given step, show how effective it is depending on which step preceded it.
    Rows are sorted by total lits removed DESC; cumulative columns show what fraction
    of the step's total work is explained by the top-N predecessor contexts."""
    if not _preproc_has_data(matched_dirs):
        return
    # Only meaningful if prev_step column exists (newer log format)
    con = sqlite3.connect("data.sqlite3")
    cur = con.cursor()
    dirs_sql = ",".join("'" + d + "'" for d in matched_dirs)
    # Check that at least some rows have a non-'none' prev_step
    cur.execute(f"SELECT COUNT(*) FROM preproc WHERE dirname IN ({dirs_sql})"
                f" AND name='{step}' AND prev_step != 'none'")
    if cur.fetchone()[0] == 0:
        con.close()
        return

    cur.execute(f"""
        SELECT prev_step,
               COUNT(*)                                                                    AS n,
               -SUM(delta_irred_long_lits)                                                 AS sum_lits,
               -SUM(delta_free_vars)                                                       AS sum_vars,
               ROUND(AVG(-delta_irred_long_lits), 1)                                       AS avg_lits,
               ROUND(AVG(-delta_free_vars), 1)                                             AS avg_vars,
               ROUND(100.0 * SUM(CASE WHEN delta_irred_long_lits < 0
                                        OR delta_free_vars < 0 THEN 1 ELSE 0 END)
                     / COUNT(*), 1)                                                        AS pct_active
        FROM preproc
        WHERE dirname IN ({dirs_sql}) AND depth = 0 AND name='{step}'
        GROUP BY prev_step
        ORDER BY sum_lits DESC
    """)
    rows = cur.fetchall()
    con.close()

    if not rows:
        return

    total_lits = sum(max(r[2] or 0, 0) for r in rows)
    total_vars = sum(max(r[3] or 0, 0) for r in rows)
    total_n    = sum(r[1] for r in rows)

    title = f"'{step}' effectiveness by predecessor step (sorted by lits removed, with cumulative)"
    print(f"\n{BLUE}{title}{RESET}")
    print(f"  Total across all predecessors: {total_n:,} invocations, "
          f"{total_lits:,} lits removed, {total_vars:,} vars removed")

    print("  (lits_rmvd/fvars_rmvd = total removed in that predecessor context;"
          " avg_lit/var = per-call average; %active = % of calls with any effect;"
          " cum_lits%/cum_vars% = cumulative share of step's total work)")
    headers = ["prev_step", "n", "%n", "lits_rmvd", "fvars_rmvd",
               "avg_lit", "avg_var", "%active", "cum_lits%", "cum_vars%"]
    cum_lits = 0
    cum_vars = 0
    str_rows = []
    for prev, n, s_lits, s_vars, a_lits, a_vars, pct_act in rows:
        s_lits = max(s_lits or 0, 0)
        s_vars = max(s_vars or 0, 0)
        cum_lits += s_lits
        cum_vars += s_vars
        pct_n      = f"{100 * n // total_n}%"       if total_n    else "0%"
        cum_lits_p = f"{100 * cum_lits // total_lits}%" if total_lits else "0%"
        cum_vars_p = f"{100 * cum_vars // total_vars}%" if total_vars else "0%"
        str_rows.append((
            str(prev),
            str(n),
            pct_n,
            f"{s_lits:,}",
            f"{s_vars:,}",
            f"{a_lits:.1f}",
            f"{a_vars:.1f}",
            f"{pct_act:.1f}%",
            cum_lits_p,
            cum_vars_p,
        ))
    _print_table(headers, str_rows)


def preproc_time_chart(matched_dirs):
    """Chart: time breakdown bar chart — total seconds per step, colored by effect type."""
    if not _preproc_has_data(matched_dirs) or not _preproc_has_step_time(matched_dirs):
        return

    dirs_sql = ",".join("'" + d + "'" for d in matched_dirs)
    con = sqlite3.connect("data.sqlite3")
    cur = con.cursor()
    cur.execute(f"""
        SELECT name,
               ROUND(SUM(step_time), 2)                                                 AS total_s,
               ROUND(SUM(CASE WHEN delta_irred_long_lits = 0 AND delta_free_vars = 0
                               AND IFNULL(delta_elimed_vars,0) = 0
                               AND IFNULL(delta_replaced_vars,0) = 0
                               AND IFNULL(delta_units,0) = 0
                          THEN step_time ELSE 0 END), 2)                                AS wasted_s
        FROM preproc
        WHERE dirname IN ({dirs_sql}) AND depth = 0 AND step_time IS NOT NULL
        GROUP BY name
        ORDER BY total_s DESC
    """)
    rows = cur.fetchall()
    con.close()

    if not rows:
        return

    lbl = _preproc_file_label(matched_dirs)
    dat_file  = f"{TMP_DIR}/preproc_time_{lbl}.dat"
    pdf_file  = f"{TMP_DIR}/preproc_time_{lbl}.pdf"
    png_file  = f"{TMP_DIR}/preproc_time_{lbl}.png"
    gp_file   = f"{TMP_DIR}/preproc_time_{lbl}.gnuplot"

    n = len(rows)
    with open(dat_file, "w") as f:
        f.write("# step  useful_s  wasted_s\n")
        for name, total_s, wasted_s in rows:
            useful_s = (total_s or 0) - (wasted_s or 0)
            f.write(f'"{name}"\t{useful_s:.2f}\t{wasted_s:.2f}\n')

    width_cm = max(22, n * 2.2)

    with open(gp_file, "w") as f:
        for term, out in [
            (f'pdfcairo size {width_cm:.0f}cm,40cm background "#d0d0d0"', pdf_file),
            (f'pngcairo size {int(width_cm * 40)},640 background "#d0d0d0"', png_file),
        ]:
            f.write(f'set terminal {term}\n')
            f.write(f'set output "{out}"\n')
            f.write(f'set title "Preprocessing time per step (useful vs wasted/no-op time)\\n{lbl}"\n')
            f.write('set ylabel "Time (seconds)"\n')
            f.write('set yrange [0:*]\n')
            f.write('set grid ytics\n')
            f.write('set key top right\n')
            f.write('set style fill solid 0.8 border -1\n')
            f.write('set style data histograms\n')
            f.write('set style histogram rowstacked\n')
            f.write('set boxwidth 0.6\n')
            f.write('set xtics rotate by -45 right\n')
            f.write('set bmargin 10\n')
            f.write(f'plot "{dat_file}" using 2:xtic(1) lc rgb "steelblue" title "useful time",\\\n')
            f.write(f'     "{dat_file}" using 3 lc rgb "red" title "no-op waste"\n\n')

    title = "Preproc time chart: useful vs no-op time per step"
    print(f"\n{BLUE}{title}{RESET}")
    print(f"  PDF: {pdf_file}  PNG: {png_file}")
    _gnuplot_run(gp_file, png_file)


def preproc_efficiency_chart(matched_dirs):
    """Chart: lits removed per second for each step (efficiency ratio)."""
    if not _preproc_has_data(matched_dirs) or not _preproc_has_step_time(matched_dirs):
        return

    dirs_sql = ",".join("'" + d + "'" for d in matched_dirs)
    con = sqlite3.connect("data.sqlite3")
    cur = con.cursor()
    cur.execute(f"""
        SELECT name,
               -SUM(delta_irred_long_lits)  AS sum_lits,
               ROUND(SUM(step_time), 2)     AS total_s
        FROM preproc
        WHERE dirname IN ({dirs_sql}) AND depth = 0 AND step_time IS NOT NULL AND step_time > 0
        GROUP BY name
        HAVING sum_lits > 0 AND total_s > 0
        ORDER BY (sum_lits / total_s) DESC
    """)
    rows = cur.fetchall()
    con.close()

    if not rows:
        return

    lbl = _preproc_file_label(matched_dirs)
    dat_file = f"{TMP_DIR}/preproc_eff_{lbl}.dat"
    pdf_file = f"{TMP_DIR}/preproc_eff_{lbl}.pdf"
    png_file = f"{TMP_DIR}/preproc_eff_{lbl}.png"
    gp_file  = f"{TMP_DIR}/preproc_eff_{lbl}.gnuplot"

    n = len(rows)
    with open(dat_file, "w") as f:
        f.write("# idx  lits_per_s  step\n")
        for i, (name, lits, total_s) in enumerate(rows):
            lits_per_s = (lits or 0) / total_s if total_s else 0
            f.write(f"{i+1}\t{lits_per_s:.0f}\t{name}\n")

    xtics = ", ".join(f'"{name}" {i+1}' for i, (name, _, __) in enumerate(rows))
    width_cm = max(20, n * 2.2)

    with open(gp_file, "w") as f:
        for term, out in [
            (f'pdfcairo size {width_cm:.0f}cm,40cm background "#d0d0d0"', pdf_file),
            (f'pngcairo size {int(width_cm * 40)},640 background "#d0d0d0"', png_file),
        ]:
            f.write(f'set terminal {term}\n')
            f.write(f'set output "{out}"\n')
            f.write(f'set title "Preprocessing efficiency: lits removed per second of CPU time\\n{lbl}"\n')
            f.write('set ylabel "Lits removed per second"\n')
            f.write('set yrange [0:*]\n')
            f.write(f'set xrange [0.5:{n + 0.5}]\n')
            f.write('set grid ytics\n')
            f.write('set key off\n')
            f.write('set style fill solid 0.8 border -1\n')
            f.write('set boxwidth 0.6\n')
            f.write(f'set xtics ({xtics}) rotate by -45 left\n')
            f.write('set bmargin 10\n')
            f.write(f'plot "{dat_file}" using ($1):2 with boxes lc rgb "steelblue" notitle\n\n')

    title = "Preproc efficiency chart: lits removed per second per step"
    print(f"\n{BLUE}{title}{RESET}")
    print(f"  PDF: {pdf_file}  PNG: {png_file}")
    _gnuplot_run(gp_file, png_file)


def preproc_time_pie_chart(matched_dirs):
    """Pie chart: share of total preprocessing wall time per step.
    Steps with <1% share are merged into 'other'."""
    if not _preproc_has_data(matched_dirs) or not _preproc_has_step_time(matched_dirs):
        return

    dirs_sql = ",".join("'" + d + "'" for d in matched_dirs)
    con = sqlite3.connect("data.sqlite3")
    cur = con.cursor()
    cur.execute(f"""
        SELECT name, ROUND(SUM(step_time), 2) AS total_s
        FROM preproc
        WHERE dirname IN ({dirs_sql}) AND depth = 0 AND step_time IS NOT NULL
        GROUP BY name
        ORDER BY total_s DESC
    """)
    rows = cur.fetchall()
    con.close()

    if not rows:
        return

    grand_total = sum(r[1] for r in rows if r[1])
    if grand_total <= 0:
        return

    # Merge steps < 1% into "other"
    kept, other_s = [], 0.0
    for name, total_s in rows:
        pct = 100.0 * (total_s or 0) / grand_total
        if pct >= 1.0:
            kept.append((name, total_s or 0))
        else:
            other_s += total_s or 0
    if other_s > 0:
        kept.append(("other", other_s))

    lbl = _preproc_file_label(matched_dirs)
    dat_file = f"{TMP_DIR}/preproc_timepie_{lbl}.dat"
    pdf_file = f"{TMP_DIR}/preproc_timepie_{lbl}.pdf"
    png_file = f"{TMP_DIR}/preproc_timepie_{lbl}.png"
    gp_file  = f"{TMP_DIR}/preproc_timepie_{lbl}.gnuplot"

    with open(dat_file, "w") as f:
        f.write("# step  seconds  pct\n")
        for name, s in kept:
            pct = 100.0 * s / grand_total
            f.write(f'"{name}"\t{s:.2f}\t{pct:.1f}\n')

    # gnuplot pie chart via filled circles — gnuplot doesn't have native pie,
    # so we use a parametric polar approach with filled boxes as wedges.
    # We write an explicit set of 'object circle' / wedge commands per slice.
    n_slices = len(kept)
    # Assign gnuplot line colors cyclically from a palette
    colors = ["#4e79a7","#f28e2b","#e15759","#76b7b2","#59a14f",
              "#edc948","#b07aa1","#ff9da7","#9c755f","#bab0ac",
              "#d37295","#fabfd2","#8cd17d","#b6992d","#499894"]

    with open(gp_file, "w") as f:
        for term, out, sz in [
            (f'pdfcairo size 20cm,16cm font ",18" background "#d0d0d0"', pdf_file, ""),
            (f'pngcairo size 800,640 font ",12" background "#d0d0d0"',   png_file, ""),
        ]:
            f.write(f'set terminal {term}\n')
            f.write(f'set output "{out}"\n')
            f.write(f'set title "Preprocessing time breakdown (total {grand_total:.0f}s)\\n{lbl}" font ",20"\n')
            f.write('set size ratio -1\n')
            f.write('unset border\n')
            f.write('unset tics\n')
            f.write('unset key\n')
            f.write('set xrange [-1.5:2.5]\n')
            f.write('set yrange [-1.3:1.3]\n')
            # Draw wedges as filled polygons
            angle = 0.0
            label_lines = []
            for i, (name, s) in enumerate(kept):
                pct = s / grand_total
                end_angle = angle + pct * 2 * math.pi
                color = colors[i % len(colors)]
                # Build polygon points for the wedge (center + arc)
                steps = max(3, int(pct * 60))
                arc_pts = []
                for k in range(steps + 1):
                    a = angle + (end_angle - angle) * k / steps
                    arc_pts.append(f"{math.cos(a):.4f},{math.sin(a):.4f}")
                poly = "0,0 to " + " to ".join(arc_pts) + " to 0,0"
                f.write(f'set object {i+1} polygon from {poly} '
                        f'fc rgb "{color}" fs solid 0.85 border lc rgb "white" lw 0.5\n')
                # Label at midpoint of arc
                mid = (angle + end_angle) / 2
                r_label = 0.65
                lx = r_label * math.cos(mid)
                ly = r_label * math.sin(mid)
                pct_str = f"{pct*100:.1f}%"
                label_lines.append(f'set label "{name}\\n{pct_str}" at {lx:.3f},{ly:.3f} center font ",13"\n')
                # Legend entry at right
                lx2 = 1.25
                ly2 = 0.95 - i * 0.13
                f.write(f'set object {n_slices+i+1} rect from {lx2-0.05},{ly2-0.04} to {lx2+0.05},{ly2+0.04} '
                        f'fc rgb "{color}" fs solid 0.85 border lc rgb "grey"\n')
                label_lines.append(f'set label "{name} ({pct*100:.1f}%)" at {lx2+0.08},{ly2} left font ",13"\n')
                angle = end_angle
            for ll in label_lines:
                f.write(ll)
            f.write('plot NaN notitle\n\n')

    title = "Preproc time pie chart: share of total preprocessing time per step"
    print(f"\n{BLUE}{title}{RESET}")
    print(f"  PDF: {pdf_file}  PNG: {png_file}")
    _gnuplot_run(gp_file, png_file)


def _display_png(png_file):
    if os.path.exists(png_file):
        with open(png_file, "rb") as fh:
            img_b64 = base64.b64encode(fh.read()).decode()
        w, h = _png_dimensions(png_file)
        print(f"\033]1337;File=inline=1;width={w}px;height={h}px:{img_b64}\a")


def _png_dimensions(png_file):
    """Read width/height from PNG header (bytes 16-24)."""
    try:
        with open(png_file, "rb") as fh:
            fh.seek(16)
            w = int.from_bytes(fh.read(4), "big")
            h = int.from_bytes(fh.read(4), "big")
        return w, h
    except Exception:
        return 800, 600


def _gnuplot_run(gp_file, png_file):
    """Run gnuplot and display the result inline."""
    os.system(f"gnuplot {gp_file}")
    _display_png(png_file)


def preproc_share_chart(matched_dirs):
    """Chart A: what fraction of total lits/vars removal each step accounts for.
    Shows two horizontal 100%-normalised bars side by side (lits and vars).
    Steps with <1% share in both are merged into 'other'."""
    if not _preproc_has_data(matched_dirs):
        return

    dirs_sql = ",".join("'" + d + "'" for d in matched_dirs)
    con = sqlite3.connect("data.sqlite3")
    cur = con.cursor()
    cur.execute(f"""
        SELECT name,
               -SUM(delta_irred_long_lits) as sum_lits,
               -SUM(delta_free_vars)       as sum_vars
        FROM preproc
        WHERE dirname IN ({dirs_sql}) AND depth = 0
        GROUP BY name
        ORDER BY sum_lits DESC
    """)
    rows = cur.fetchall()
    con.close()

    total_lits = sum(max(r[1], 0) for r in rows)
    total_vars = sum(max(r[2], 0) for r in rows)
    if total_lits == 0 and total_vars == 0:
        return

    # Keep steps with >= 1% share in either metric; merge rest into "other"
    kept, other_lits, other_vars = [], 0, 0
    for name, lits, vars_ in rows:
        lits = max(lits, 0); vars_ = max(vars_, 0)
        pct_l = 100.0 * lits / total_lits if total_lits else 0
        pct_v = 100.0 * vars_ / total_vars if total_vars else 0
        if pct_l >= 1.0 or pct_v >= 1.0:
            kept.append((name, lits, vars_))
        else:
            other_lits += lits; other_vars += vars_
    if other_lits > 0 or other_vars > 0:
        kept.append(("other", other_lits, other_vars))

    # Sort kept by lits% descending so largest bar is on the left
    kept.sort(key=lambda r: r[1], reverse=True)
    n = len(kept)
    width_cm = max(20, n * 2.0)

    lbl = _preproc_file_label(matched_dirs)
    dat_file = f"{TMP_DIR}/preproc_share_{lbl}.dat"
    pdf_file = f"{TMP_DIR}/preproc_share_{lbl}.pdf"
    png_file = f"{TMP_DIR}/preproc_share_{lbl}.png"
    gp_file  = f"{TMP_DIR}/preproc_share_{lbl}.gnuplot"

    with open(dat_file, "w") as f:
        f.write("# idx  pct_lits  pct_vars  step_name\n")
        for i, (name, lits, vars_) in enumerate(kept):
            pl = 100.0 * lits / total_lits if total_lits else 0
            pv = 100.0 * vars_ / total_vars if total_vars else 0
            f.write(f"{i+1}\t{pl:.1f}\t{pv:.1f}\t{name}\n")

    xtics = ", ".join(f'"{name}" {i+1}' for i, (name, _, __) in enumerate(kept))
    with open(gp_file, "w") as f:
        for term, out in [
            (f'pdfcairo size {width_cm:.0f}cm,40cm background "#d0d0d0"', pdf_file),
            (f'pngcairo size {int(width_cm * 40)},640 background "#d0d0d0"', png_file),
        ]:
            f.write(f'set terminal {term}\n')
            f.write(f'set output "{out}"\n')
            f.write(f'set title "Share of total preprocessing work per step\\n{lbl}"\n')
            f.write('set ylabel "% of all lits / vars removed across all CNFs"\n')
            f.write('set yrange [0:100]\n')
            f.write(f'set xrange [0.5:{n + 0.5}]\n')
            f.write('set grid ytics\n')
            f.write('set key top right\n')
            f.write('set style fill solid 0.7 border -1\n')
            f.write('set boxwidth 0.3\n')
            f.write(f'set xtics ({xtics}) rotate by -45 left\n')
            f.write('set bmargin 10\n')
            f.write(f'plot "{dat_file}" using ($1-0.18):2 with boxes lc rgb "steelblue" title "% lits removed",\\\n')
            f.write(f'     "{dat_file}" using ($1+0.18):3 with boxes lc rgb "dark-green" title "% vars removed"\n\n')

    title = "Preproc share chart: % of total lits/vars removed per step"
    print(f"\n{BLUE}{title}{RESET}")
    print(f"  PDF: {pdf_file}  PNG: {png_file}")
    _gnuplot_run(gp_file, png_file)



def preproc_cumulative_chart(matched_dirs):
    """Chart C: cumulative SUM of lits/vars removed across pipeline stages (ordered by avg position).
    Each step is one horizontal tick; lines show how the running total grows."""
    if not _preproc_has_data(matched_dirs):
        return

    dirs_sql = ",".join("'" + d + "'" for d in matched_dirs)
    con = sqlite3.connect("data.sqlite3")
    cur = con.cursor()

    # Per-step: sum removed and average pipeline position (avg step_num)
    cur.execute(f"""
        SELECT name,
               -SUM(delta_irred_long_lits) as sum_lits,
               -SUM(delta_irred_long_cls)  as sum_cls,
               -SUM(delta_irred_bins)      as sum_bins,
               -SUM(delta_free_vars)       as sum_vars,
               SUM(IFNULL(delta_units, 0)) as sum_units,
               AVG(step_num)               as avg_pos
        FROM preproc
        WHERE dirname IN ({dirs_sql}) AND depth = 0
        GROUP BY name
        ORDER BY avg_pos ASC
    """)
    rows = cur.fetchall()
    con.close()

    if not rows:
        return

    lbl = _preproc_file_label(matched_dirs)
    dat_file = f"{TMP_DIR}/preproc_cumul_{lbl}.dat"
    pdf_file = f"{TMP_DIR}/preproc_cumul_{lbl}.pdf"
    png_file = f"{TMP_DIR}/preproc_cumul_{lbl}.png"
    gp_file  = f"{TMP_DIR}/preproc_cumul_{lbl}.gnuplot"

    # Normalise to millions for readability
    cum_lits = cum_cls = cum_vars = cum_units = 0.0
    with open(dat_file, "w") as f:
        f.write("# i  cum_lits_M  cum_cls_M  cum_vars_M  cum_units_M  step_name\n")
        f.write("0\t0.0\t0.0\t0.0\t0.0\tstart\n")
        for i, (name, s_lits, s_long_cls, s_bins, s_vars, s_units, _pos) in enumerate(rows):
            cum_lits  += max(s_lits, 0) / 1e6
            cum_cls   += max((s_long_cls or 0) + (s_bins or 0), 0) / 1e6
            cum_vars  += max(s_vars, 0) / 1e6
            cum_units += max(s_units or 0, 0) / 1e6
            f.write(f"{i+1}\t{cum_lits:.3f}\t{cum_cls:.3f}\t{cum_vars:.3f}\t{cum_units:.3f}\t{name}\n")

    n = len(rows)
    height_cm = max(8, n * 0.40)
    width_cm  = 24

    ytics_parts = ['"start" 0'] + [f'"{name}" {i+1}' for i, (name, *_) in enumerate(rows)]
    ytics_str = ", ".join(ytics_parts)

    with open(gp_file, "w") as f:
        for term, out in [
            (f'pdfcairo size {width_cm}cm,{height_cm:.0f}cm background "#d0d0d0"', pdf_file),
            (f'pngcairo size {width_cm * 50},{int(height_cm * 50)} background "#d0d0d0"', png_file),
        ]:
            f.write(f'set terminal {term}\n')
            f.write(f'set output "{out}"\n')
            f.write(f'set title "Cumulative preprocessing effect across pipeline (total across all CNFs)\\n{lbl}"\n')
            f.write('set xlabel "Cumulative total removed (millions)"\n')
            f.write('set ylabel ""\n')
            f.write(f'set yrange [-0.5:{n + 0.5}]\n')
            f.write('set xrange [0:*]\n')
            f.write('set grid xtics\n')
            f.write('set key bottom right\n')
            f.write(f'set ytics ({ytics_str})\n')
            f.write(f'plot "{dat_file}" using 2:1 with linespoints lc rgb "steelblue" lw 2 pt 7 ps 1 title "lits removed",\\\n')
            f.write(f'     "{dat_file}" using 3:1 with linespoints lc rgb "dark-orange" lw 2 pt 5 ps 1 title "cls removed (bin+long)",\\\n')
            f.write(f'     "{dat_file}" using 4:1 with linespoints lc rgb "dark-green" lw 2 pt 9 ps 1 title "vars removed",\\\n')
            f.write(f'     "{dat_file}" using 5:1 with linespoints lc rgb "dark-violet" lw 2 pt 11 ps 1 title "units found"\n\n')

    title = f"Preproc cumulative chart (all steps, n={n}, ordered by pipeline position)"
    print(f"\n{BLUE}{title}{RESET}")
    print(f"  PDF: {pdf_file}  PNG: {png_file}")
    _gnuplot_run(gp_file, png_file)


def print_preproc_per_step_detail(matched_dirs, verbose=False):
    """Top steps by total absolute lits impact (all occurrences), sorted by median lits reduction."""
    if not _preproc_has_data(matched_dirs):
        return

    dirs_sql = ",".join("'" + d + "'" for d in matched_dirs)
    con = sqlite3.connect("data.sqlite3")
    stats = _preproc_step_stats(con, dirs_sql)

    # Also need n_active_any (active on any metric) — fetch from DB
    cur = con.cursor()
    cur.execute(
        f"SELECT name,"
        f" ROUND(100.0 * SUM(CASE WHEN delta_irred_long_lits != 0 OR delta_irred_bins != 0"
        f"   OR delta_free_vars != 0 THEN 1 ELSE 0 END) / COUNT(*), 1)"
        f" FROM preproc WHERE dirname IN ({dirs_sql}) AND depth = 0 GROUP BY name"
    )
    pct_active = {row[0]: row[1] for row in cur.fetchall()}
    con.close()

    if not stats:
        return

    # Sort by median lits reduction (most reduction first)
    ordered = sorted(stats.items(), key=lambda kv: kv[1]['med_lits'])[:10]

    title = "Top steps by median lits reduction (all occurrences)"
    print(f"\n{BLUE}{title}{RESET}")

    headers = ["step", "n", "med_d_lits", "med_d_bins", "med_d_vars", "%active"]
    str_rows = [
        (name,
         str(s['n']),
         str(s['med_lits']),
         str(s['med_bins']),
         str(s['med_vars']),
         str(pct_active.get(name, 0)))
        for name, s in ordered
    ]
    _print_table(headers, str_rows)


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
        dat_file = f"{TMP_DIR}/scatter_{safe1}_vs_{safe2}.dat"
        pdf_file = f"{TMP_DIR}/scatter_{safe1}_vs_{safe2}.pdf"
        png_file = f"{TMP_DIR}/scatter_{safe1}_vs_{safe2}.png"
        gp_file  = f"{TMP_DIR}/scatter_{safe1}_vs_{safe2}.gnuplot"

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
                (f'pdfcairo size 15cm,15cm background "#d0d0d0"', pdf_file),
                (f'pngcairo size 600,600 background "#d0d0d0"',   png_file),
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
        _display_png(png_file)

    con.close()


def generate_gnuplot(fname2_s, verbose=False):
    gnuplotfn = f"{TMP_DIR}/cdf.gnuplot"
    pdf_file = f"{TMP_DIR}/cdf.pdf"
    png_file = f"{TMP_DIR}/cdf.png"
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
            ('pdfcairo size 45cm,65cm background "#d0d0d0"', pdf_file),
            ('pngcairo size 600,600 background "#d0d0d0"',   png_file),
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

    filename = f'{TMP_DIR}/overview.ipynb'
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
    # "out-ganak-mccomp2324-1229753-0", # lots of bug fixes, beauty changes with Claude, etc
    # "out-ganak-mccomp2324-1231407-0", # the same as above but without (most) of the Claude improvements
    # "out-ganak-mccomp2324-1247484-0", ## new shrinking, fixing arjun SLOW_DEBUG, improved propagation idea from CaDiCaL, improved 3-tier clause database. Also undoing only 2 particular Claude changes (propagation and --prob 0 non-zeroing of data)
    # "out-ganak-mccomp2324-1250247-", # CMS cleanup, oracle improvements, fix parsing issue (?) of CNF header -- weird + one of them is a binary with: reason-side bumping and lbd update, evsids
    # "out-ganak-mccomp2324-1256426-0", # fixing oracle, mostly, and also header parsing more lax
    # "out-ganak-mccomp2324-1261017-0", # Different order in Arjun
    # "out-ganak-mccomp2324-1279568-0", # release
    # "out-ganak-mccomp2324-1281478-0", # more data
    # "out-ganak-mccomp2324-1282000-0", # new preproc VERY GOOD !!!!!!!!!!!!!!!!!
    # "out-ganak-mccomp2324-1282412-0", # new preproc v2
    # "out-ganak-mccomp2324-1286085-", # 4 new puura orders
    # "out-ganak-mccomp2324-1286556-", # 4 new puura orders, and wl-based cache compact
    # "out-ganak-mccomp2324-1294423-0", # improve CMS's subsumption/strengthening, fix oracle to work non-backbone, add units from oracle, allow to set linit to cadiback
    # "out-ganak-mccomp2324-1296621-", # same as above, but fixing CMS bug
    # "out-ganak-mccomp2324-1299016-0", # trying to get back to old ganak speed -- CMS changes (reverts) GOOOOOOOD!!!!
    #"out-ganak-mccomp2324-1315264-", # Vivification, new order for backw, new setup for backw -- GOTO BUG!
    # "out-ganak-mccomp2324-1318712-3", # fixed goto1 and now vivif is turned on/off with backwtype=0
    # "out-ganak-mccomp2324-1351541-", # let's fix it, plus llm improve it(?) -- NOOOOO, we left the CMS nonsense in!
    # "out-ganak-mccomp2324-1358071-1", # same as above, but with CMS changes reverted -- Current best (kinda...)
    # "out-ganak-mccomp2324-1373930-1", # no grow after
    # "out-ganak-mccomp2324-1358071-", # same as above, but with CMS changes reverted
    # "out-ganak-mccomp2324-1452294-1", # same as above, but with arjun changes back to "out-ganak-mccomp2324-1299016-0" (with option for 1 new idea)
    #"out-ganak-mccomp2324-1458467-", # let's try some new ideas from LLM for arjun improvement -- but cadiback was wrongly set up
    # "out-ganak-mccomp2324-1514564-", # let's try some new ideas from LLM for arjun improvement -- there's been another f*ck up
    # "out-ganak-mccomp2324-1517017-0", # 4 full runs for all the 4 new arjun orders
    # "out-ganak-mccomp2324-1635700-0", # fix the printing of the preproc data
    # "out-ganak-mccomp2324-1743408", # ddnnf
    # "out-ganak-mccomp2324-1747186-0", # faster ddnnf, new hash function

    # 5 min timeout runs:
    # "out-ganak-mccomp2324-1755057-0", # 5 min timeout
    # "out-ganak-mccomp2324-1755057-3", # 5 min timeout
    # "out-ganak-mccomp2324-1758343-5", # new 5 min timeout run
    # "out-ganak-mccomp2324-1762059-", # new 5 min timeout run
    # "out-ganak-mccomp2324-1783926-2", # also extend
    # best is: --fast --tditers 100--arjunextendmaxconfl 3000
    # 0b4881b4_11e203ea_67c5648a_5e1ee18e

    # final MCC
    "out-ganak-mccomp2324-1783906-0", # final competition stuff: norm and trying kitten. Slowdown is purely machine failure/CPU overload
                                      # running ganak_0b4881b4_11e203ea_67c5648a_5e1ee18e
    # "out-ganak-mccomp2324-1812040-0", # 2 min timeout
    # "out-ganak-mccomp2324-1812431-4", # 2 min timeout, more configs
    # "out-ganak-mccomp2324-1812683-", # 2 min timeout, more configs
    ## other stuff
    # "out-ganak-mccomp2324-1783906-1", # kitten
    # "out-ganak-mccomp2324-1817408-0", # gates-eq + replace in the middle after gates-based eq
    "out-ganak-mccomp2324-1835807-", # --rdbclstarget check, running ganak_0b4881b4_11e203ea_67c5648a_5e1ee18e
]
# only_dirs = [
#      "mei-march-2026-1239767-1", # gpmc
#      # "mei-march-2026-1239767-0", # ganak old
#      # "mei-march-2026-1269673-0", # ganak release, but WRONG SED
#      "mei-march-2026-1274973-0", # ganak release, new SED
# ]

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
    parser.add_argument("--nopreproc", action="store_true",
                        help="Skip all preprocessing tables and graphs (preproc table)")
    parser.add_argument("--nopairwise", action="store_true",
                        help="No pairwise comparisons")
    parser.add_argument("--nodistribution", action="store_true",
                        help="Don't print distributions of metrics")
    parser.add_argument("--cdf", action="store_true",
                        help="ONLY generate the PAR2/solved summary table and the CDF graph; skip everything else")
    args = parser.parse_args()

    os.makedirs(TMP_DIR, exist_ok=True)

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
    if not args.cdf and not args.nopairwise:
      scatter_plot_time_pairs(matched_dirs, fname_like, args.verbose)
    fname2_s, table_todo = build_csv_data(todo, matched_dirs, only_calls, not_calls, not_versions, fname_like, args.verbose)

    if args.cdf:
        if args.verbose:
            print("Printing summary tables...")
        print_summary_tables(table_todo, fname_like, args.full, args.verbose)
        if args.verbose:
            print("Generating gnuplot script...")
        gnuplotfn, pdf_file, png_file = generate_gnuplot(fname2_s, args.verbose)
        for path in [pdf_file, png_file]:
            if os.path.exists(path):
                os.unlink(path)
        os.system(f"gnuplot {gnuplotfn}")
        console_title = "CDF: instances counted vs. solve time"
        print(f"\n{BLUE}{console_title}{RESET}")
        print(f"  PDF: {pdf_file}  PNG: {png_file}")
        _display_png(png_file)
        return

    if args.verbose:
        print(f"Selected {len(table_todo)} dir/version combinations")
    seen = set()
    for dir, _ in table_todo:
        if dir not in seen:
            seen.add(dir)
            os.system(f"./cache_miss_bucket_summary.py {dir}")

    if not args.nodistribution:
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
    if not args.nopreproc:
        print_preproc_diffs(table_todo, fname_like, args.verbose)

        # Preprocessing step analysis — one block per directory
        for d in matched_dirs:
            one = [d]
            print(f"\n{GREEN}{'='*70}{RESET}")
            print(f"{GREEN}  Preprocessing analysis  ---  {d}{RESET}")
            print(f"{GREEN}{'='*70}{RESET}")
            print_preproc_delta_table(one, args.verbose)
            print_preproc_step_efficiency(one)

            print_preproc_time_breakdown(one)
            print_preproc_noop_waste(one)
            print_preproc_extended_vars(one)
            print_step_predecessor_effectiveness(one, step="must-scc-vrepl")
            preproc_cumulative_chart(one)
            preproc_time_pie_chart(one)

    if not args.nopairwise:
      unique_dirs = list(dict.fromkeys(d for d, _ in table_todo))
      for dir1, dir2 in itertools.combinations(unique_dirs, 2):
          print_two_dir_diffs(dir1, dir2, fname_like, args.verbose)
          print_solved_only_diffs(dir1, dir2, fname_like, args.verbose)
          print_solution_mismatches(dir1, dir2, fname_like, args.verbose)

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
    _display_png(png_file)


if __name__ == "__main__":
    main()
