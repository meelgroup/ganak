from pathlib import Path
from collections import Counter
import csv,json,os,sqlite3
from datetime import datetime,timezone
P=Path(os.environ['GANAK_OUTPUT'])
def write(name,fields,rows):
    out=P/name;t=out.with_suffix('.tmp')
    with t.open('w') as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
    t.replace(out)
def publish():
    con=sqlite3.connect(P/'study.sqlite');con.row_factory=sqlite3.Row
    data=con.execute('select * from jobs order by id').fetchall();con.close()
    screen=[];pairs=[];times=[];unique=set();counts=Counter(r['state'] for r in data)
    for r in data:
        if not r['stock'] and r['state']!='error':continue
        arms={a:json.loads(r[a]) for a in ('stock','sweep','combined') if r[a]}
        s=arms.get('stock',{});screen.append(dict(id=r['id'],source=r['source'],state=r['state'],core=r['core'],stock_status=s.get('status',''),stock_s=s.get('wall',''),detail=r['detail'] or s.get('diagnostic','')))
        if r['state']!='done':continue
        assert r['sig'] and r['sig'] not in unique;unique.add(r['sig'])
        assert 60<=s['wall']<=600 and len(arms)==3 and s['count'].isdigit()
        for v in arms.values():
            assert v['status'] in ('timeout','solved') and (v['status']!='solved' or v['count']==s['count'])
        p=dict(id=r['id'],source=r['source'],family=r['family'],core=r['core'],sha256=r['sha'],token_sha256=r['sig'],count=s['count'])
        for a,v in arms.items():p[a+'_s']=v['wall'];p[a+'_status']=v['status']
        pairs.append(p)
    cov=[]
    for cap in (60,120,180,300,600):
        sets={a:{p['id'] for p in pairs if p[a+'_status']=='solved' and p[a+'_s']<=cap} for a in ('stock','sweep','combined')}
        cov.append(dict(cap_s=cap,n=len(pairs),**{a:len(v) for a,v in sets.items()},gained_vs_sweep=' '.join(sorted(sets['combined']-sets['sweep'])),lost_vs_sweep=' '.join(sorted(sets['sweep']-sets['combined']))))
    write('screen.csv',['id','source','state','core','stock_status','stock_s','detail'],screen)
    write('pairs.csv',['id','source','family','core','sha256','token_sha256','count','stock_s','stock_status','sweep_s','sweep_status','combined_s','combined_status'],pairs)
    write('coverage.csv',cov[0].keys(),cov)
    inband=sum(x['stock_status']=='solved' and 60<=x['stock_s']<=600 for x in screen)
    stamp=datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')
    lines=[f'# Ganak medium-duration benchmark',f'{stamp}. {len(pairs)} unique qualifying formulas have complete, count-checked comparisons toward the target of 200. {len(screen)} screening records are available, including {inband} in-band baseline runs before deduplication and completion of both comparisons.',
           'Stock Ganak develop 9a05d35 defines eligibility. Existing sweep bd4cb7e and the frozen combined policy use the same pinned physical core as stock for each formula. The timeout is 600 wall seconds per method and the address-space limit is 9 GB per process. Shorter budgets are computed from these same runs. A candidate timeout stays in the denominator. Exact answers must agree with stock.',
           'The local pool contains 1,000 competition IDs from 2022–2026 and 2,505 entries from the 2024 submission archive. In-band inputs are deduplicated by the DIMACS token stream, ignoring comments and whitespace. Sampling prioritizes historical difficulty and families yielding qualifying inputs. Treat coverage as conditional on the selected formulas. Core types vary and the host is shared.',
           '| Seconds | Stock | Existing sweep | Combined |','|---:|---:|---:|---:|',*[f"| {r['cap_s']} | {r['stock']} | {r['sweep']} | {r['combined']} |" for r in cov],
           'Screening states: '+', '.join(f'{k} {v}' for k,v in sorted(counts.items()))+'.']
    (P/'STATUS.tmp').write_text('\n\n'.join(lines[:4])+'\n\n'+'\n'.join(lines[4:])+'\n');(P/'STATUS.tmp').replace(P/'STATUS.md')
    return counts
if __name__=='__main__':
    import sys,time
    if '--monitor' in sys.argv:
        while True:
            try:
                counts=publish()
            except OSError as e:
                if e.errno not in (28,122):raise
                time.sleep(30)
                continue
            if not counts['running'] and (counts['done']>=200 or not counts['pending']):break
            time.sleep(30)
    else:print(publish())
