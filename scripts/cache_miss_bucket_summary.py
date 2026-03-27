#!/usr/bin/env python3
"""Summarise solver metrics by cache-miss-rate bucket for a given dirname.

Prints a table with columns for high/medium/low cache miss rate, showing
median values for key metrics. Run from the directory containing mydb.sql.
"""

import argparse
import sqlite3


def median(cur, col, where):
    cur.execute(f"SELECT COUNT(*) FROM data WHERE {where} AND {col} IS NOT NULL")
    n = cur.fetchone()[0]
    if n == 0:
        return None
    cur.execute(
        f"SELECT {col} FROM data WHERE {where} AND {col} IS NOT NULL"
        f" ORDER BY {col} LIMIT 1 OFFSET {n // 2}"
    )
    return cur.fetchone()[0]


def print_table(rows, headers):
    widths = [max(len(headers[i]), max(len(r[i]) for r in rows)) for i in range(len(headers))]
    sep = "+-" + "-+-".join("-" * w for w in widths) + "-+"
    fmt = "| " + " | ".join(f"{{:<{w}}}" for w in widths) + " |"
    print(sep)
    print(fmt.format(*headers))
    print(sep)
    for row in rows:
        print(fmt.format(*row))
    print(sep)


def main():
    parser = argparse.ArgumentParser(description="Cache miss rate bucket summary")
    parser.add_argument("dirname", help="dirname to analyse (must exist in mydb.sql)")
    parser.add_argument("--cutoff", type=float, default=100.0,
                        help="Minimum ganak_time - arjun_time to include (default: 100)")
    args = parser.parse_args()

    con = sqlite3.connect("mydb.sql")
    cur = con.cursor()

    cur.execute("SELECT COUNT(*) FROM data WHERE dirname=?", (args.dirname,))
    if cur.fetchone()[0] == 0:
        cur.execute("SELECT DISTINCT dirname FROM data ORDER BY dirname")
        known = [r[0] for r in cur.fetchall()]
        con.close()
        parser.error(f"dirname '{args.dirname}' not found in mydb.sql. Known dirnames:\n  "
                     + "\n  ".join(known))

    base = (f"dirname='{args.dirname}'"
            f" AND cache_miss_rate IS NOT NULL"
            f" AND (ganak_time - arjun_time) >= {args.cutoff}")

    buckets = [
        ("high (>=0.85)",     f"{base} AND cache_miss_rate >= 0.85"),
        ("medium (0.4-0.85)", f"{base} AND cache_miss_rate >= 0.4 AND cache_miss_rate < 0.85"),
        ("low (<0.4)",        f"{base} AND cache_miss_rate < 0.4"),
    ]

    metrics = [
        ("num instances",       None),
        ("median td_width",     "td_width"),
        ("median compsK",       "compsK"),
        ("median indep_sz",     "indep_sz"),
        ("median mem_mb",       "ganak_mem_mb"),
        ("median solve_t (s)",  "ganak_time - arjun_time"),
        ("median conflicts M",  "conflicts / 1000000.0"),
    ]

    rows = []
    for metric, col in metrics:
        row = [metric]
        for _, where in buckets:
            if col is None:
                cur.execute(f"SELECT COUNT(*) FROM data WHERE {where}")
                val = str(cur.fetchone()[0])
            else:
                val = median(cur, col, where)
                if val is None:
                    val = "N/A"
                elif isinstance(val, float):
                    val = f"{val:.2f}"
                else:
                    val = f"{val:,}"
            row.append(val)
        rows.append(row)

    cur.execute("SELECT count(*), sum(ganak_time IS NOT NULL) FROM data WHERE dirname=?", (args.dirname,))
    total_rows, total_solved = cur.fetchone()
    con.close()

    headers = ["metric"] + [b[0] for b in buckets]
    print(f"\nCache miss rate bucket summary for: {args.dirname}")
    print(f"NOTE: filtered to instances where ganak_time - arjun_time >= {args.cutoff}s")
    print(f"      (i.e. hard instances only; {total_solved} solved out of {total_rows} total are NOT all shown)")
    print()
    print_table(rows, headers)


if __name__ == "__main__":
    main()
