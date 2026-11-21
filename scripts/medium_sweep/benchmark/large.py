"""Resumable local screen. Only stock runtime and input metadata influence selection."""
from bench import P,CORES,runarm,cnf_stream
import collections,contextlib,csv,fcntl,gzip,hashlib,json,lzma,math,os,re,resource,shutil,sqlite3,subprocess,sys,tarfile,tempfile,time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,as_completed
archive=Path(os.environ['GANAK_DATA'])/'heur/reference/collected.tar'
def db():
    c=sqlite3.connect(P/'study.sqlite',timeout=120,isolation_level=None)
    c.row_factory=sqlite3.Row
    return c
def freeze():
    c=db();c.execute('pragma journal_mode=WAL')
    c.execute('create table if not exists jobs(id text primary key,source text,family text,size int,priority int,ord text,state text default "pending",core int,sha text,sig text,stock text,sweep text,combined text,detail text)')
    if c.execute('select count(*) from jobs').fetchone()[0]:return
    with gzip.open(Path(__file__).with_name('manifest.csv.gz'), 'rt') as f:
        c.execute('begin')
        c.executemany('insert into jobs(id,source,family,size,priority,ord) values(:id,:source,:family,:size,:priority,:ord)', csv.DictReader(f))
    c.execute('commit')
@contextlib.contextmanager
def stream(r):
    if r['source']=='track':
        y,i=map(int,r['id'].split('-'))
        with cnf_stream(y,i) as f:yield f
    else:
        with tarfile.open(archive) as t,t.extractfile(r['source']) as packed:
            if r['source'].endswith('.xz'):
                with lzma.open(packed) as f:yield f
            elif r['source'].endswith('.gz'):
                with gzip.open(packed) as f:yield f
            else:yield packed
def signature(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for line in f:
            words=line.split()
            if words and words[0]!=b'c':h.update(b' '.join(words)+b' ')
    return h.hexdigest()
def claim(c,cpu):
    c.execute('begin immediate')
    if c.execute("select count(*) from jobs where state='done'").fetchone()[0]>=200:
        c.execute('commit');return
    # Prioritize any interrupted pair on its original core.
    r=c.execute("select * from jobs where state='pending' and core=? order by ord limit 1",(cpu,)).fetchone()
    if not r:
        rows=c.execute("select * from jobs where state='pending' and core is null order by priority,ord").fetchall()
        if rows:
            ceiling=rows[0]['priority'];pool=[x for x in rows if x['priority']==ceiling]
            completed=c.execute("select family,count(*) as trials,sum(case when state='done' then 1 else 0 end) as wins from jobs where state in ('fast','done','timeout','duplicate','error') group by family").fetchall()
            stats={x['family']:(x['trials'],x['wins']) for x in completed}
            issued=c.execute("select count(*) from jobs where state!='pending'").fetchone()[0]
            def score(x):
                n,k=stats.get(x['family'],(0,0))
                return (k+1)/(n+6)
            r=pool[0] if issued%4==0 else max(pool,key=score)
    if r:c.execute("update jobs set state='running',core=? where id=?",(cpu,r['id']))
    c.execute('commit');return r
def lane(cpu):
    os.sched_setaffinity(0,{cpu});resource.setrlimit(resource.RLIMIT_AS,(9_000_000_000,9_000_000_000));c=db()
    while True:
        avail=int(next(x.split()[1] for x in Path('/proc/meminfo').read_text().splitlines() if x.startswith('MemAvailable:')))
        if avail<6*1024**2 or shutil.disk_usage(P).free<8<<20:time.sleep(5);continue
        r=claim(c,cpu)
        if not r:return
        n=r['id']
        try:
            with os.fdopen(os.memfd_create('scale200-input.cnf'), 'wb') as backing:
                path=Path(f'/proc/{os.getpid()}/fd/{backing.fileno()}');h=hashlib.sha256()
                with stream(r) as f,path.open('wb') as out:
                    while chunk:=f.read(1<<20):
                        h.update(chunk);out.write(chunk)
                s=json.loads(r['stock']) if r['stock'] else runarm('stock',path,600,n)
                c.execute('update jobs set stock=?,sha=? where id=?',(json.dumps(s),h.hexdigest(),n))
                if s['status']!='solved' or s['wall']<60:
                    state='fast' if s['status']=='solved' else s['status']
                    c.execute('update jobs set state=? where id=?',(state,n));continue
                sig=signature(path);c.execute('begin immediate')
                dup=c.execute("select id from jobs where sig=? and id!=? and state not in ('error','duplicate') limit 1",(sig,n)).fetchone()
                c.execute('update jobs set sig=? where id=?',(sig,n))
                if dup:c.execute("update jobs set state='duplicate',detail=? where id=?",(dup['id'],n))
                c.execute('commit')
                if dup:continue
                okay=True
                for a in (('sweep','combined') if int(r['ord'][0],16)%2 else ('combined','sweep')):
                    v=json.loads(r[a]) if r[a] else runarm(a,path,600,n)
                    c.execute('update jobs set '+a+'=? where id=?',(json.dumps(v),n))
                    if v['status']=='error' or (v['status']=='solved' and v['count']!=s['count']):okay=False
                c.execute('update jobs set state=? where id=?',('done' if okay else 'error',n))
        except Exception as e:
            if c.in_transaction:c.execute('rollback')
            c.execute("update jobs set state='error',detail=? where id=?",(repr(e),n))
if __name__=='__main__':
    lock=(P/'run.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    freeze();db().execute("update jobs set state='pending' where state='running'")
    with ProcessPoolExecutor(len(CORES)) as ex:
        for f in as_completed([ex.submit(lane,cpu) for cpu in CORES]):f.result()
