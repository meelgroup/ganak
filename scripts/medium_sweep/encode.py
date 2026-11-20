from itertools import combinations

def packet(n, clauses):
    owners={};groups=[];pairs=set()
    for c in clauses:
        if 2<=len(c)<=63 and min(c)>0:
            if len(groups)>=128:return
            v=len(groups);groups.append(c)
            for x,b in enumerate(c):
                if b in owners or b>n:return
                owners[b]=(v,x)
        elif len(c)==2 and max(c)<0 and c[0]!=c[1]:
            pairs.add(tuple(sorted((-c[0],-c[1]))))
        else:return
    if len(owners)!=n or not n:return
    if any((min(a,b),max(a,b)) not in pairs for g in groups for a,b in combinations(g,2)):return
    cross=[]
    for a,b in sorted(pairs):
        if a not in owners or b not in owners:return
        v,x=owners[a];u,y=owners[b]
        if v!=u:cross.append(f'{v} {x} {u} {y}\n')
    return (f'{len(groups)} {len(cross)}\n'+' '.join(str(len(g)) for g in groups)+'\n'+''.join(cross)).encode()
