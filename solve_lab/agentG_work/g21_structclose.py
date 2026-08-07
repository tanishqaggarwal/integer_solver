"""SOUND closure: structural (modp=None) free-input support, so DEAD monomials
(second-order activation paths) are not missed."""
import os, sys, json, time
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentG_work')
import gsym, gGclose
from gsym import *
import suppfree
SRC=sys.argv[1] if len(sys.argv)>1 else '/home/user/integer_solver/solve_lab/s10/AG_39013.json'
FL=[int(x) for x in sys.argv[2].split(',') if x] if len(sys.argv)>2 and sys.argv[2]!='-' else []
OUT=sys.argv[3] if len(sys.argv)>3 else 'closed_struct.json'
v=L.load(SRC); ad.fwd(v,rounds=6)
for b in FL: v[b]=1-v[b]
ad.fwd(v,rounds=8)
vm=[x%P for x in v]
idx,freelist,vs = suppfree.build(vm, modp=None)
FREESET=set(freelist)
def supp(a):
    m=0
    for w in L.avars[a]: m |= vs[w] if w<len(vs) else 0
    s={freelist[i] for i in range(len(freelist)) if (m>>i)&1}
    s |= {w for w in L.avars[a] if w in FREESET}
    return s
av=L.all_atom_values(v)
S=set()
for a in gsym.check_atoms():
    if av[a]%P: S |= {u for u in supp(a) if not gGclose.isbool(u)}
print('seed |S| =',len(S), flush=True)
for it in range(20):
    SY=sorted(S); n=len(SY)
    t0=time.time()
    try:
        val=gsym.build(v,SY,cap=6,verbose=False)
    except OverflowError:
        print('  degree cap hit at |S|=%d'%n); break
    sym=[]
    maxd=0; tot=0
    for a in gsym.check_atoms():
        f=gsym.evalpoly_sym(a,val,n,6)
        if not isinstance(f,int) or f%P:
            sym.append(a)
            if not isinstance(f,int): maxd=max(maxd,gsym.deg(f)); tot+=len(f)
    new=set()
    for a in sym: new |= {u for u in supp(a) if not gGclose.isbool(u)}
    print('  it%d |S|=%d sym=%d maxdeg=%d terms=%d new=%d (%.1fs)'%(it,n,len(sym),maxd,tot,len(new-S),time.time()-t0), flush=True)
    if new<=S: break
    S|=new
json.dump(sorted(S), open(OUT,'w'))
print('STRUCTURAL CLOSED SET (%d) -> %s'%(len(S),OUT))
