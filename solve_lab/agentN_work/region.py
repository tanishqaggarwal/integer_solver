"""Residual region in the agent-H frame: which atoms live ENTIRELY inside it (compensators)."""
import ev, fast, json, os, time
from fast import St, csup, inv
from close2 import close
from collections import defaultdict
p=115792089237316195423570985008687907853269984665640564039457584007908834671663
BITS=json.load(open('bits.json')); ALL=set(BITS['A']+BITS['B'])
eq_of=defaultdict(list)
for i,(m,sq,tl) in enumerate(ev.eq_terms):
    for c,a in tl: eq_of[a].append(i)
st0=St({})
b=BITS['A'][0]
st=st0.clone().set_free({b:1})
out,ok,_=close(st,frozen=set(ALL),maxsteps=300)
print('state: score',out.score(),'nz',sorted(out.nz()),'failing',len(out.fails))
NZ=sorted(out.nz())
E=set()
for a in NZ: E.update(eq_of[a])
print('region |E| =',len(E),'  currently failing',len(out.fails),' satisfied in region',len(E)-len(out.fails))
# atoms whose ENTIRE equation footprint lies inside E
comp=[a for a in ev.checks if a not in NZ and eq_of[a] and set(eq_of[a])<=E]
print('atoms fully inside the region (candidate compensators):',len(comp),comp[:40])
# which of them are settable without touching anything outside E?
AV=ev.atom_vars; FREE=set(ev.F['free0'])
settable=[]
for a in comp+NZ:
    for X in AV[a]:
        if X not in FREE: continue
        g=out.clone().set_free({X:out.fv.get(X,0)+1})
        moved = g.av[a]!=out.av[a]
        outside = any(g.av[q]!=out.av[q] for q in g.av if q not in set(comp)|set(NZ))
        if moved and not outside:
            settable.append((a,X)); break
print('settable inside-region atoms:',len(settable),settable[:30])
json.dump({'NZ':NZ,'E':sorted(E),'comp':comp,'settable':[[a,x] for a,x in settable]},open('region.json','w'))
