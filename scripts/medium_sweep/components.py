from collections import Counter

def factors(n,clauses):
    parent=list(range(n+1))
    def find(x):
        while parent[x]!=x:
            parent[x]=parent[parent[x]];x=parent[x]
        return x
    for c in clauses:
        representative=find(abs(c[0]))
        for x in c:parent[find(abs(x))]=representative
    parts={}
    for c in clauses:parts.setdefault(find(abs(c[0])),[]).append(c)
    keys=Counter()
    for part in parts.values():
        ids={v:i for i,v in enumerate(sorted({abs(x) for c in part for x in c}),1)}
        key=tuple(sorted({tuple(sorted(ids[x] if x>0 else -ids[-x] for x in c)) for c in part}))
        keys[key]+=1
    return keys
