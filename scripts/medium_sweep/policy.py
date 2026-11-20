"""Count an ordinary model-counting DIMACS file with a bounded structural prepass."""
import argparse
import re
import subprocess
import os
import tempfile
import time
from pathlib import Path
from algo import load, simplify
from components import factors
from equivalent import reduce
P=Path(__file__).resolve().parent
from encode import packet


def count(path,cap,ganak=None):
    ganak=ganak or os.environ.get("GANAK", "ganak")
    start=time.monotonic()
    def left():return max(0,cap-(time.monotonic()-start))
    n,clauses=load(path);free=0
    with tempfile.TemporaryDirectory() as directory:
        if any(len(c)==1 for c in clauses):
            result=reduce(n,clauses)
            if result is None:return 0
            v,cs,removed,_,_=result
            separated=False
            if v>.75*n:
                # Keep the original clause order after unit propagation if it
                # exposes two sizable independent components.
                unc=simplify(n,clauses)
                assert unc is not None
                uv,uc,uf,_=unc
                parts=factors(uv,uc)
                large=sum(mult for key,mult in parts.items()
                          if len({abs(x) for c in key for x in c})*20>=uv)
                if large>=2:
                    v,cs,removed=uv,uc,uf
                    separated=True
            if separated or v<=.75*n:
                free=removed
                if v==0:return 1<<free
                data=''.join(' '.join(map(str,c))+' 0\n' for c in cs)
                if all((len(c)==2 and min(c)>0) or (len(c)>2 and max(c)<0) for c in cs):
                    try:
                        p=subprocess.run([P/'enumerate'],input=f'{v} {len(cs)}\n'+data,text=True,capture_output=True,
                                         timeout=min(5,left()),check=True)
                        if p.stdout.strip().isdigit():return int(p.stdout)<<free
                    except subprocess.TimeoutExpired:pass
                path=Path(directory)/'reduced.cnf'
                path.write_text(f'p cnf {v} {len(cs)}\n'+data)
        elif n<=2000 and all((len(c)==2 and max(c)<0) or (c and min(c)>0) for c in clauses):
            z=packet(n, clauses)
            if z:
                p=subprocess.run([P/'resume','adaptive',str(1<<24)],input=z,
                                 capture_output=True,check=True,timeout=left())
                answer=p.stdout.split()
                if answer[0]==b'exact':return int(answer[1])
        p=subprocess.run([ganak,'--threads','1','--prob','0','--frontierstates','8192','--groupedsweep','1',
                          '--frontierminvars','4','--frontierprofile','1',path],capture_output=True,timeout=left())
        match=re.search(rb'^c s exact arb int (\d+)\s*$',p.stdout,re.M)
        if p.returncode not in (0,10,20) or not match:raise RuntimeError('Ganak failed to return an exact count')
        return int(match[1])<<free


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('cnf',type=Path)
    parser.add_argument('--ganak',help='Sweep-enabled Ganak executable (defaults to $GANAK or ganak)')
    parser.add_argument('--seconds',type=float,default=600);args=parser.parse_args()
    if args.seconds<=0:parser.error('Use a positive timeout')
    start=time.monotonic()
    try:
        result=count(args.cnf,args.seconds,args.ganak)
        print('c s exact arb int',result) if time.monotonic()-start<=args.seconds else print('s UNKNOWN')
    except subprocess.TimeoutExpired:print('s UNKNOWN')
