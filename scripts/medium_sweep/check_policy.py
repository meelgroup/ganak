from policy import count
from algo import simplify
from pathlib import Path
from tempfile import TemporaryDirectory
import random
def oracle(n,cs):
    return sum(all(any(bool(x>0)==bool(mask&(1<<(abs(x)-1))) for x in c) for c in cs) for mask in range(1<<n))

rng=random.Random(688115)
with TemporaryDirectory() as d:
    p=Path(d)/'input.cnf'
    for i in range(100):
        n=rng.randint(3,9)
        if i%4==0:
            cs=[(-1,2),(1,-2),(3,),(-3,1)]
        elif i%4==1:
            n=8;groups=[(x,x+1) for x in (1,3,5,7)]
            cs=groups+[(-a,-b) for a,b in groups]+[(-a,-b) for j,g in enumerate(groups) for h in groups[j+1:] for a in g for b in h if rng.random()<.18]
        else:
            cs=[tuple(rng.choice((-1,1))*rng.randint(1,n) for _ in range(rng.randint(1,3))) for _ in range(rng.randint(1,16))]
        p.write_text(f'p cnf {n} {len(cs)}\n'+''.join(' '.join(map(str,c))+' 0\n' for c in cs))
        v=oracle(n,cs)
        assert count(p,15)==v
        z=simplify(n,cs,True);assert (0 if z is None else oracle(z[0],z[1])<<z[2])==v
print('100 whole-policy runs and 100 literal-order-preserving reductions matched direct enumeration.')
