import os, sys, time
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentG_work')
import gsym2 as G
from gsym2 import L, ad, P
import gGclose
SRC=os.environ.get('SRC','/home/user/integer_solver/solve_lab/s10/AG_39013.json')
FL=[int(x) for x in sys.argv[1].split(',') if x] if sys.argv[1]!='-' else []
TARG=[int(x) for x in sys.argv[2].split(',')]
v=L.load(SRC); ad.fwd(v,rounds=6)
for b in FL: v[b]=1-v[b]
ad.fwd(v,rounds=8)
FREE=[u for u in range(L.NVARS) if u not in L.definer]
NB=[u for u in FREE if not gGclose.isbool(u)]
val,sk=G.build(v,NB,cap=6)
def mstr(m): return '*'.join('x%d'%NB[k]+('^%d'%e if e>1 else '') for k,e in m) or '1'
def sc(c):
    c%=P
    return str(c) if c<10**12 else ('(-%d)'%(P-c) if P-c<10**12 else 'C[%s..]'%str(c)[:10])
print('frame',FL,'score',L.NEQ-len(L.failing_eqs(L.all_atom_values(v))))
for a in TARG:
    f=G.evalatom(a,val,6)
    if isinstance(f,int):
        print('a%-6d (%d eqs): CONSTANT %s'%(a,len(L.atom2eq.get(a,{})),sc(f))); continue
    print('a%-6d (%d eqs) deg%d terms%d:'%(a,len(L.atom2eq.get(a,{})),G.deg(f),len(f)))
    for m,c in sorted(f.items(),key=lambda kv:(-G.mdeg(kv[0]),kv[0]))[:25]:
        print('     %-40s %s'%(mstr(m),sc(c)))
