import sqlite3
import os
os.chdir("/home/soos/development/sat_solvers/ganak/build/data")
con = sqlite3.connect("data.sqlite3")
cur = con.cursor()

cur.execute("""
SELECT name,
       COUNT(DISTINCT fname) as n_cnfs,
       SUM(CASE WHEN delta_irred_long_lits < 0 THEN 1 ELSE 0 END) as invoc_lits,
       SUM(CASE WHEN delta_free_vars < 0 THEN 1 ELSE 0 END) as invoc_vars,
       -SUM(delta_irred_long_lits) as sum_lits,
       -SUM(delta_free_vars) as sum_vars,
       COUNT(*) as n_invoc
FROM preproc
GROUP BY name
ORDER BY sum_lits DESC
""")
rows = cur.fetchall()
print(f"{'step':<30} {'n_cnfs':>7} {'%inv_lits':>9} {'%inv_vars':>9} {'sum_lits':>12} {'sum_vars':>12} {'n_invoc':>8}")
print("-"*104)
for r in rows:
    name, n_cnfs, invoc_lits, invoc_vars, sum_lits, sum_vars, n_invoc = r
    pct_lits = 100*invoc_lits/n_invoc if n_invoc else 0
    pct_vars = 100*invoc_vars/n_invoc if n_invoc else 0
    print(f"{name:<30} {n_cnfs:>7} {pct_lits:>8.0f}% {pct_vars:>8.0f}% {sum_lits:>12,} {sum_vars:>12,} {n_invoc:>8,}")

cur.execute("SELECT -SUM(delta_irred_long_lits), -SUM(delta_free_vars) FROM preproc")
t = cur.fetchone()
print(f"\nTOTAL: lits={t[0]:,}  vars={t[1]:,}")

print("\n\n--- Per-CNF distribution of lits removed (steps that remove lits) ---")
cur.execute("""
SELECT p.fname,
       -SUM(p.delta_irred_long_lits) as total_lits_removed,
       -SUM(p.delta_free_vars) as total_vars_removed,
       d.ganak_time
FROM preproc p
JOIN data d ON d.fname = p.fname AND d.dirname = p.dirname
WHERE d.ganak_time IS NOT NULL
GROUP BY p.dirname, p.fname
ORDER BY total_lits_removed DESC
LIMIT 20
""")
rows2 = cur.fetchall()
print(f"{'fname':<50} {'lits_removed':>13} {'vars_removed':>13} {'time':>8}")
for r in rows2:
    print(f"{r[0][-48:]:<50} {r[1]:>13,} {r[2]:>13,} {r[3]:>8.2f}")

print("\n\n--- Steps that are ALWAYS active (invoc_pct_lits > 80%) vs rarely active ---")
cur.execute("""
SELECT name,
       COUNT(*) as n_invoc,
       -SUM(delta_irred_long_lits) as sum_lits,
       ROUND(100.0 * SUM(CASE WHEN delta_irred_long_lits < 0 THEN 1 ELSE 0 END) / COUNT(*), 1) as pct_lits
FROM preproc
GROUP BY name
HAVING sum_lits > 0
ORDER BY pct_lits DESC
""")
rows3 = cur.fetchall()
print(f"{'step':<30} {'n_invoc':>8} {'sum_lits':>12} {'%active':>8}")
for r in rows3:
    print(f"{r[0]:<30} {r[1]:>8,} {r[2]:>12,} {r[3]:>7.1f}%")

print("\n\n--- occ-bve deep dive: bins added ---")
cur.execute("""
SELECT -SUM(delta_irred_bins), SUM(delta_irred_bins), COUNT(*)
FROM preproc WHERE name = 'occ-bve'
""")
r = cur.fetchone()
print(f"delta_bins total={r[0]:,} sum={r[1]:,} n={r[2]:,}")

print("\n\n--- What's removing the most lits ON AVERAGE per invocation? ---")
cur.execute("""
SELECT name,
       COUNT(*) as n_invoc,
       -SUM(delta_irred_long_lits) as sum_lits,
       ROUND(-1.0*SUM(delta_irred_long_lits)/COUNT(*), 1) as avg_lits_per_invoc,
       ROUND(-1.0*SUM(delta_free_vars)/COUNT(*), 1) as avg_vars_per_invoc
FROM preproc
GROUP BY name
ORDER BY avg_lits_per_invoc DESC
""")
rows4 = cur.fetchall()
print(f"{'step':<30} {'n_invoc':>8} {'sum_lits':>12} {'avg_lits/inv':>13} {'avg_vars/inv':>13}")
for r in rows4:
    print(f"{r[0]:<30} {r[1]:>8,} {r[2]:>12,} {r[3]:>13,.1f} {r[4]:>13,.1f}")
