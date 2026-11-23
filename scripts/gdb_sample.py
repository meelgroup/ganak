# Poor man's profiler for when perf record is not allowed (perf_event_paranoid >= 1).
# usage: PMP_ARGS="--mode 1 f.cnf" PMP_N=40 PMP_DT=0.5 PMP_DELAY=10 PMP_OUT=s.txt gdb -batch -x scripts/gdb_sample.py build/ganak
#        then: python3 scripts/gdb_sample.py --agg s.txt
import os, sys, re, collections, subprocess

def aggregate(fname):
    samples = open(fname).read().split('=== sample')[1:]
    top = collections.Counter(); incl = collections.Counter()
    for s in samples:
        fr = re.findall(r'^#\d+\s+(?:0x[0-9a-f]+ in )?([^\s(]+)', s, re.M)
        if not fr: continue
        top[fr[0]] += 1
        for f in set(fr): incl[f] += 1
    print(len(samples), "samples")
    print("self:"); [print(f"  {v:4d} {k}") for k, v in top.most_common(20)]
    print("inclusive:"); [print(f"  {v:4d} {k}") for k, v in incl.most_common(30)]

if len(sys.argv) > 2 and sys.argv[1] == '--agg':
    aggregate(sys.argv[2])
    sys.exit(0)

import gdb
args = os.environ['PMP_ARGS']
n = int(os.environ.get('PMP_N', '40'))
dt = float(os.environ.get('PMP_DT', '1'))
delay = float(os.environ.get('PMP_DELAY', '0'))
outname = os.environ['PMP_OUT']
out = open(outname, 'w')
gdb.execute("set pagination off")
gdb.execute("set debuginfod enabled off")
gdb.execute("set print frame-arguments none")
# noprint would imply nostop
gdb.execute("handle SIGUSR1 stop print nopass")
gdb.execute("starti " + args + " > " + outname + ".prog 2>&1")
pid = gdb.selected_inferior().pid
for i in range(n):
    subprocess.Popen(["sh", "-c", f"sleep {delay if i == 0 else dt}; kill -USR1 {pid}"])
    try: gdb.execute("continue")
    except gdb.error: break
    if not gdb.selected_inferior().pid: break
    out.write("=== sample %d\n%s" % (i, gdb.execute("bt 25", to_string=True)))
    out.flush()
gdb.execute("kill")
