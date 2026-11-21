from pathlib import Path
import gzip,os,re,shutil,signal,subprocess,sys,time
from archives import cnf_stream
P=Path(os.environ['GANAK_OUTPUT']).resolve()
P.mkdir(exist_ok=True,parents=True)
CORES=list(map(int,os.environ['GANAK_CORES'].split(',')))
bins={'stock':os.environ['GANAK_STOCK'],'sweep':os.environ['GANAK']}
def runarm(name,path,cap,case_id):
    cmd=[bins['stock' if name=='stock' else 'sweep'],'--threads','1','--prob','0']
    if name=='sweep':cmd+=['--frontierstates','8192','--groupedsweep','1','--frontierminvars','4','--frontierprofile','1']
    if name=='combined':cmd=[sys.executable,str(Path(__file__).resolve().parents[1]/'policy.py'),'--seconds',str(cap)]
    cmd.append(str(path))
    start=time.monotonic()
    with os.fdopen(os.memfd_create('scale200-output'), 'w+b') as log:
        pr=subprocess.Popen(cmd,stdout=log,stderr=log,start_new_session=True)
        timed_out=False
        try:pr.wait(timeout=cap)
        except subprocess.TimeoutExpired:
            timed_out=True
            try:os.killpg(pr.pid,signal.SIGKILL)
            except ProcessLookupError:pass
            pr.wait()
        wall=time.monotonic()-start
        log.seek(0)
        folder=P/'logs'
        folder.mkdir(exist_ok=True)
        target=folder/f'{case_id}-{name}.log.gz'
        tmp=target.with_suffix('.tmp')
        with gzip.open(tmp,'wb') as out:
            shutil.copyfileobj(log,out)
        tmp.replace(target)
        if timed_out or wall>cap:
            return dict(status='timeout',wall=wall,count='')
        log.seek(0)
        count=None
        unknown=False
        for line in log:
            if line.startswith(b'c s exact arb int '):
                m=re.fullmatch(rb'c s exact arb int (\d+)\s*',line)
                if m:count=m[1].decode()
            if line.strip()==b's UNKNOWN':unknown=True
        if count is not None and pr.returncode in (0,10,20):
            return dict(status='solved',wall=wall,count=count)
        if unknown:return dict(status='timeout',wall=wall,count='')
        log.seek(max(0,log.tell()-1600))
        return dict(status='error',wall=wall,count='',diagnostic=log.read().decode(errors='replace'))
