"""Eliminate signed equivalences proved by binary implications."""
from algo import simplify

def collapse(n,clauses):
    adj=[[] for _ in range(2*n)];rev=[[] for _ in adj]
    def key(x):return 2*(abs(x)-1)+(x<0)
    for clause in clauses:
        if len(clause)!=2:continue
        a,b=map(key,clause)
        for u,v in ((a^1,b),(b^1,a)):
            adj[u].append(v);rev[v].append(u)
    seen=bytearray(2*n);order=[]
    for root in range(2*n):
        if seen[root]:continue
        stack=[root]
        while stack:
            x=stack.pop()
            if x<0:order.append(~x);continue
            if seen[x]:continue
            seen[x]=1;stack.append(~x)
            stack.extend(v for v in adj[x] if not seen[v])
    comp=[-1]*(2*n)
    for u in reversed(order):
        if comp[u]>=0:continue
        stack=[u];comp[u]=u
        while stack:
            for v in rev[stack.pop()]:
                if comp[v]<0:comp[v]=u;stack.append(v)
    labels={};mapped=[0]
    for i in range(n):
        a,b=comp[2*i:2*i+2]
        if a==b:return None
        if a not in labels:
            label=1+len(labels)//2;labels[a]=label;labels[b]=-label
        mapped.append(labels[a])
    size=len(labels)//2
    if size==n:return n,clauses
    return size,[tuple(mapped[x] if x>0 else -mapped[-x] for x in c) for c in clauses]

def reduce_units(initial):
    if initial is None:
        return None
    n, clauses, free, units = initial
    merged = 0
    while True:
        result = collapse(n, clauses)
        if result is None:
            return None
        size, clauses = result
        if size == n:
            return n, clauses, free, units, merged
        merged += n - size
        result = simplify(size, clauses)
        if result is None:
            return None
        n, clauses, f, u = result
        free += f
        units += u
