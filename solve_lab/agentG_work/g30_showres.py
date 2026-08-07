import os, sys, pickle
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentG_work')
import g29_frame as F
import gsym2 as G
from gsym2 import L, ad, P
arg=sys.argv[1]
FL=[int(x) for x in arg.split(',') if x] if arg!='-' else []
v=list(F.v0)
for b in FL: v[b]=1-v[b]
ad.fwd(v,rounds=8)
r=F.analyse(v)
NB=F.NB
def mstr(m): return '*'.join('x%d'%NB[k]+('^%d'%e if e>1 else '') for k,e in m) or '1'
def sc(c):
    c%=P
    return str(c) if c<10**12 else ('(-%d)'%(P-c) if P-c<10**12 else 'C[%s..]'%str(c)[:10])
print('flip',FL,'rank',r['rank'],'nfree',r['nfree'],'ninc',len(r['incchecks']),'nzc',len(r['nzc']))
vars_=set()
for a,g in r['res']:
    print('\na%-6d (%d eqs):'%(a,len(L.atom2eq.get(a,{}))))
    if isinstance(g,int): print('   CONST',g%P); continue
    for m,c in sorted(g.items(),key=lambda kv:(-G.mdeg(kv[0]),kv[0])):
        print('   %-40s %s'%(mstr(m),sc(c)))
        for k,e in m: vars_.add(NB[k])
print('\nvariables in the residual system:',sorted(vars_))
pickle.dump({'res':r['res'],'NB':NB,'flip':FL},open('res_%s.pkl'%(arg.replace(',','_')),'wb'))
