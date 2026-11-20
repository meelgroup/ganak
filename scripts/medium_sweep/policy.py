"""Bounded structural prepass for ordinary unweighted DIMACS."""
import argparse
import gzip
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import threading
import time

from algo import load, simplify
from components import factors
from encode import packet
from equivalent import reduce_units

P = Path(__file__).resolve().parent
PROFILE = (b'c frontier ', b'c grouped ', b'c sweep_profile ', b'c o Ganak SHA1:')


def ganak_count(cmd, left, trace):
    # Forward diagnostics as they arrive, so an external kill preserves progress.
    # Only the outer driver knows how many free variables to restore.
    answer = []
    def read_output(pipe):
        for line in pipe:
            if line.startswith(b'c s exact arb int '):
                match = re.fullmatch(rb'c s exact arb int (\d+)\s*', line)
                if match:
                    answer.append(match[1])
            elif trace and line.startswith(PROFILE):
                print(line.decode(errors='replace').rstrip(), flush=True)

    with subprocess.Popen(cmd, stdout=subprocess.PIPE) as child:
        reader = threading.Thread(target=read_output, args=(child.stdout,))
        reader.start()
        try:
            child.wait(timeout=left())
        finally:
            if child.poll() is None:
                child.kill()
            child.wait()
            reader.join()
        if child.returncode not in (0, 10, 20) or len(answer) != 1:
            raise RuntimeError('Ganak failed to return a single exact count')
    return int(answer[0])


def count(path, cap, ganak=None, *, trace=False, buckets=True, previous=False, kernel_timing=False, maxcache=None):
    start = time.monotonic()
    def left():
        return max(0, cap - (time.monotonic() - start))
    def event(name, **fields):
        if trace:
            values = ' '.join(f'{key}={value}' for key, value in fields.items())
            print(f'c medium event={name} wall_s={time.monotonic()-start:.6f} {values}', flush=True)

    event('begin', policy=2, reduction='previous' if previous else 'conservative',
          buckets=int(buckets), kernel_timing=int(kernel_timing), node_cap=1 << 24)
    with tempfile.TemporaryDirectory() as directory:
        with path.open('rb') as raw:
            compressed = raw.read(2) == b'\x1f\x8b'
        if compressed:
            with gzip.open(path, 'rb') as src:
                path = Path(directory)/'input.cnf'
                with path.open('wb') as dest:
                    shutil.copyfileobj(src, dest)
        n, clauses = load(path)
        if left() <= 0:
            raise subprocess.TimeoutExpired('prepass', cap)
        event('parsed', variables=n, clauses=len(clauses))
        free = 0
        route = 'original'
        if any(len(c) == 1 for c in clauses):
            units = simplify(n, clauses)
            reduced = reduce_units(units)
            if reduced is None:
                event('finish', route='conflict')
                return 0
            v, cs, removed, fixed, merged = reduced
            uv, uc, uf, uu = units
            event('reduced', unit_vars=uv, reduced_vars=v, fixed=fixed, merged=merged, free=removed)
            use_units = not previous and uv <= .75*n and v > uv/2
            if v > .75*n:
                parts = factors(uv, uc)
                large = sum(mult for key, mult in parts.items()
                            if len({abs(x) for c in key for x in c})*20 >= uv)
                use_units = large >= 2
            if use_units or v <= .75*n:
                if use_units:
                    v, cs, removed = uv, uc, uf
                route = 'units' if use_units else 'reduced'
                free = removed
                if v == 0:
                    event('finish', route='fully_reduced', free=free)
                    return 1 << free
                data = ''.join(' '.join(map(str,c))+' 0\n' for c in cs)
                if all((len(c)==2 and min(c)>0) or (len(c)>2 and max(c)<0) for c in cs):
                    event('probe', method='enumerate')
                    try:
                        p = subprocess.run([P/'enumerate'], input=f'{v} {len(cs)}\n'+data,
                                           text=True, capture_output=True, timeout=min(5,left()), check=True)
                        if p.stdout.strip().isdigit():
                            event('finish', route='enumerate', free=free)
                            return int(p.stdout) << free
                        event('decline', method='enumerate', reason='model_cap')
                    except subprocess.TimeoutExpired:
                        event('decline', method='enumerate', reason='timeout')
                path = Path(directory)/'reduced.cnf'
                path.write_text(f'p cnf {v} {len(cs)}\n'+data)
        elif n<=2000 and all((len(c)==2 and max(c)<0) or (c and min(c)>0) for c in clauses):
            z = packet(n, clauses)
            if z:
                event('probe', method='adaptive_domain')
                p = subprocess.run([P/'resume','adaptive',str(1<<24)], input=z,
                                   capture_output=True, check=True, timeout=left())
                status, value, nodes, *_ = p.stdout.split()
                if status == b'exact':
                    event('finish', route='adaptive_domain', nodes=int(nodes))
                    return int(value)
                event('decline', method='adaptive_domain', reason='node_cap', nodes=int(nodes))

        event('ganak', route=route, free=free)
        cmd = [ganak or os.environ.get('GANAK','ganak'),'--threads','1','--prob','0',
               '--frontierstates','8192','--groupedsweep','1','--frontierminvars','4','--frontierprofile','2' if kernel_timing else '1']
        if maxcache is not None:
            cmd += ['--maxcache',str(maxcache)]
        if not buckets:
            cmd += ['--frontierbucket','0']
        value = ganak_count(cmd+[path], left, trace)
        event('finish', route=route, free=free)
        return value << free


if __name__=='__main__':
    # Competition counts can exceed Python's default decimal conversion limit.
    if hasattr(sys, 'set_int_max_str_digits'):
        sys.set_int_max_str_digits(0)
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('cnf', type=Path)
    parser.add_argument('--ganak', help='Sweep-enabled Ganak executable (defaults to $GANAK or ganak)')
    parser.add_argument('--seconds', type=float, default=600)
    parser.add_argument('--kernel-timing', action='store_true')
    parser.add_argument('--maxcache', type=int, help='Forward the native cache limit to Ganak')
    parser.add_argument('--no-buckets', action='store_true')
    parser.add_argument('--previous-reduction', action='store_true')
    args = parser.parse_args()
    if args.seconds <= 0:
        parser.error('Use a positive timeout')
    start = time.monotonic()
    try:
        result = count(args.cnf, args.seconds, args.ganak, trace=True, buckets=not args.no_buckets,
                       previous=args.previous_reduction, kernel_timing=args.kernel_timing, maxcache=args.maxcache)
        if time.monotonic()-start > args.seconds:
            raise subprocess.TimeoutExpired('policy', args.seconds)
        print('c s exact arb int', result)
    except subprocess.TimeoutExpired:
        print('c medium event=timeout', flush=True)
        print('s UNKNOWN')
