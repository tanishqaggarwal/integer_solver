import os, sys, json
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentG_work')
import gsym
from gsym import *
import suppfree
src = '/home/user/integer_solver/solve_lab/s10/AG_39013.json'
v = L.load(src); ad.fwd(v, rounds=6); vm=[x%P for x in v]
TARG=[688,1618,19297,19299,30984,36185,40608,40812]
idx, freelist, vs = suppfree.build(vm, modp=True)
tot=set()
for a in TARG:
    m=0
    for w in L.avars[a]: m |= vs[w] if w<len(vs) else 0
    tot |= {freelist[i] for i in range(len(freelist)) if (m>>i)&1}
    tot |= {w for w in L.avars[a] if w not in L.definer}
# boolean detection: exists atom whose poly is c*(u*u) - c*(u)  (any scaling)
BOOL=set()
for u in tot:
    for a in L.var_atoms[u]:
        pl=L.polys[a]
        if len(pl)==2:
            ks=list(pl.keys())
            if sorted(map(len,ks))==[1,2]:
                m1=[k for k in ks if len(k)==1][0]; m2=[k for k in ks if len(k)==2][0]
                if m1==(u,) and m2==(u,u) and pl[m1]==-pl[m2]:
                    BOOL.add(u); break
NB=sorted(tot-BOOL)
print('union %d ; boolean %d ; NON-BOOLEAN %d' % (len(tot),len(BOOL),len(NB)))
print('non-boolean:', NB)
for u in NB: print('   x%-6d bits=%d  val%%p=%s' % (u, v[u].bit_length(), (v[u]%P) if v[u].bit_length()<40 else str(v[u]%P)[:20]+'...'))
json.dump({'union':sorted(tot),'bool':sorted(BOOL),'nonbool':NB}, open('supp8.json','w'))
