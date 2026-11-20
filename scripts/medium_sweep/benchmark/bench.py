from pathlib import Path
import os,re,signal,subprocess,sys,time
from archives import cnf_stream
P=Path(os.environ['GANAK_OUTPUT']).resolve()
P.mkdir(exist_ok=True,parents=True)
CORES=list(map(int,os.environ['GANAK_CORES'].split(',')))
bins={'stock':os.environ['GANAK_STOCK'],'sweep':os.environ['GANAK']}
def runarm(name,path,cap):
    cmd=[bins['stock' if name=='stock' else 'sweep'],'--threads','1','--prob','0']
    if name=='sweep':cmd+=['--frontierstates','8192','--groupedsweep','1','--frontierminvars','4','--frontierprofile','1']
    if name=='combined':cmd=[sys.executable,str(Path(__file__).resolve().parents[1]/'policy.py'),'--seconds',str(cap)]
    cmd.append(str(path))
    start=time.monotonic()
    with os.fdopen(os.memfd_create('scale200-output'), 'w+b') as log:
        pr=subprocess.Popen(cmd,stdout=log,stderr=log,start_new_session=True)
        try:pr.wait(timeout=cap)
        except subprocess.TimeoutExpired:
            try:os.killpg(pr.pid,signal.SIGKILL)
            except ProcessLookupError:pass
            pr.wait()
            return dict(status='timeout',wall=time.monotonic()-start,count='')
        wall=time.monotonic()-start;log.seek(0,2);log.seek(max(0,log.tell()-65536));out=log.read()
    m=re.search(rb'^c s exact arb int (\d+)\s*$',out,re.M)
    if m and pr.returncode in (0,10,20) and wall<=cap:return dict(status='solved',wall=wall,count=m[1].decode())
    if wall>=cap or b's UNKNOWN' in out:return dict(status='timeout',wall=wall,count='')
    return dict(status='error',wall=wall,count='',diagnostic=out[-1600:].decode(errors='replace'))
