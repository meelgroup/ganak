"""Exact preprocessing for plain, unweighted, unprojected DIMACS."""
from collections import defaultdict, deque

def load(path):
    clauses=[];pending=[];n=m=None
    for line in path.read_text().splitlines():
        tok=line.split()
        if not tok:continue
        if tok[0]=='c':
            if len(tok)>1 and tok[1].lower() in ('ind','p','show','weight','vp'):raise ValueError('annotated input')
            if tok[1:2]==['t'] and tok[2:]!=['mc']:raise ValueError('unsupported track')
            continue
        if tok[:2]==['p','cnf']:n,m=map(int,tok[2:]);continue
        for t in map(int,tok):
            if t:pending.append(t)
            else:clauses.append(tuple(pending));pending=[]
    assert not pending and len(clauses)==m and all(0<abs(x)<=n for c in clauses for x in c)
    return n,clauses

def simplify(n,clauses,keep_order=False):
    occurrences=defaultdict(list);queue=deque();rows=[];used=set();original=[]
    for c in clauses:
        s=set(c)
        if any(-x in s for x in s):continue
        if not s:return None
        if len(s)==1:queue.append(next(iter(s)))
        for x in s:occurrences[x].append(len(rows))
        rows.append(s)
        if keep_order:original.append(c)
    while queue:
        x=queue.popleft()
        if -x in used:return None
        if x in used:continue
        used.add(x)
        for i in occurrences[x]:rows[i]=None
        for i in occurrences[-x]:
            s=rows[i]
            if s is None:continue
            s.discard(-x)
            if not s:return None
            if len(s)==1:queue.append(next(iter(s)))
    cs=[s for s in rows if s];active=sorted({abs(x) for s in cs for x in s});ids={v:i+1 for i,v in enumerate(active)}
    residual=[tuple(sorted(ids[abs(x)] if x>0 else -ids[-x] for x in s)) for s in cs]
    if keep_order:
        residual=[tuple(ids[x] if x>0 else -ids[-x] for x in dict.fromkeys(c) if x in s) for c,s in zip(original,rows) if s]
    return len(active),residual,n-len(used)-len(active),len(used)
