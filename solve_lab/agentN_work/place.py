"""Where can the two-congruence defect be placed?  Measure every combination."""
import ev, fast, json
from fast import St
from close2 import close
from collections import defaultdict
p=115792089237316195423570985008687907853269984665640564039457584007908834671663
C1=125787314747601108116039725163361763116550465675981151838811516827327919228823597744635626
C2=91416258160755509149180373473728639746431157665678710450404458852172057265575180278101002
BITS=json.load(open('bits.json')); ALL=set(BITS['A']+BITS['B'])
eq_of=defaultdict(list)
for i,(m,sq,tl) in enumerate(ev.eq_terms):
    for c,a in tl: eq_of[a].append(i)
st0=St({}); b=BITS['A'][0]
base=st0.clone().set_free({b:1})
out,_,_=close(base,frozen=set(ALL),maxsteps=300)
V37892=out.v[37892]; V13682=out.v[13682]
print('delivered x_37892 mod p == C1 ?',V37892%p==C1%p)
print('delivered x_13682 mod p == C2 ?',V13682%p==C2%p)
res=[]
for n1,v1 in (('a30982 clean (x_18956:=x_37892)',V37892),('a688 clean (x_18956:=C1)',C1)):
    for n2,v2 in (('a30980 clean (x_24468:=x_13682)',V13682),('a1618 clean (x_24468:=C2)',C2)):
        st=st0.clone().set_free({b:1,18956:v1,24468:v2})
        o,_,_=close(st,frozen=set(ALL)|{18956,24468},maxsteps=300)
        E=set()
        for a in o.nz(): E.update(eq_of[a])
        res.append((o.score(),len(o.fails),len(E),sorted(o.nz()),n1,n2))
        print('%-38s | %-36s -> score %d  failing %d  region %d  nz %s'%(n1,n2,o.score(),len(o.fails),len(E),sorted(o.nz())))
res.sort(reverse=True)
print('\nBEST placement in this frame: score',res[0][0])
