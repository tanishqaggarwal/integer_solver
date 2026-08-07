import os, sys, json
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentG_work')
import gsym, gred, gGclose
from gsym import *
import suppfree
SRC='/home/user/integer_solver/solve_lab/s10/AG_39013.json'
FL=[int(x) for x in sys.argv[1].split(',')] if len(sys.argv)>1 else [2081,24601]
v0=L.load(SRC); ad.fwd(v0,rounds=6)
v=list(v0)
for b in FL: v[b]=1-v[b]
ad.fwd(v,rounds=8)
S=gGclose.closure(v)
r=gred.reduce_state(v,S)
print('flip',FL,'|S|',len(S),'rank',r['rank'],'ninc',r['ninc'],'nzc',len(r['nzc']))
vm=[x%P for x in v]
idx,freelist,vs=suppfree.build(vm,modp=True)
FREESET=set(freelist)
def supp(a):
    m=0
    for w in L.avars[a]: m |= vs[w] if w<len(vs) else 0
    s={freelist[i] for i in range(len(freelist)) if (m>>i)&1}
    s |= {w for w in L.avars[a] if w in FREESET}
    return s
for a,val in r['nzc']:
    sp=supp(a)
    nb=[u for u in sp if u not in gGclose._BC or not gGclose._BC[u]]
    print('\na%-6d neq=%d val=%d'%(a,len(L.atom2eq.get(a,{})),val))
    print('   supp size %d ; supp=%s'%(len(sp),sorted(sp)[:40]))
    for m,c in sorted(L.polys[a].items(), key=lambda kv:(-len(kv[0]),kv[0]))[:12]:
        vs2=' * '.join('x%d[%s]'%(u,str(v[u])[:10]) for u in m)
        print('     %-60s coeff %s'%(vs2, c if abs(c)<10**12 else str(c)[:12]+'..'))
