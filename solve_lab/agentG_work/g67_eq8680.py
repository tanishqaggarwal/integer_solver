"""Locate eq8680 in my own equation-level model: which unknowns it involves, its degree,
its value at x0, and whether it pins an unknown to zero.  Also eq29125."""
import os, sys, pickle, collections
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentG_work')
import gsym2 as G
from gsym2 import L, ad, P
D=pickle.load(open('/home/user/integer_solver/solve_lab/agentG_work/coset_model.pkl','rb'))
Lin=pickle.load(open('/home/user/integer_solver/solve_lab/agentG_work/coset_lin.pkl','rb'))
NB=D['NB']; lin=D['lin']; non=D['non']; x0=Lin['x0']; pt=D['pt']; col=Lin['col']
n=len(NB); ix={u:i for i,u in enumerate(NB)}
allq={i:('LIN',f) for i,f in lin}
for i,f in non: allq[i]=('NON',f)
def ev(f,vec):
    if isinstance(f,int): return f%P
    s=0
    for m,c in f.items():
        t=c
        for k,e in m: t=t*pow(vec[k],e,P)%P
        s=(s+t)%P
    return s
for q in [8680,29125,12231,12270,12350,14584,18673,22044,2554,9421]:
    if q not in allq:
        print('eq%-6d : identically zero in my model (no unknown reaches it)'%q); continue
    kind,f=allq[q]
    if isinstance(f,int):
        print('eq%-6d : %s constant %d'%(q,kind,f%P)); continue
    vs=sorted({NB[k] for m in f for k,_ in m})
    print('eq%-6d : %s deg %d, %d terms, %d unknowns %s'%(q,kind,G.deg(f),len(f),len(vs),vs[:12]))
    print('          value at x0 = %s ; at the deliverable = %s'
          %('0' if ev(f,x0)==0 else 'NONZERO','0' if ev(f,pt)==0 else 'NONZERO'))
    if G.deg(f)<=2 and len(vs)<=6:
        for m,c in sorted(f.items(),key=lambda kv:-G.mdeg(kv[0])):
            nm='*'.join('x%d'%NB[k]+('^%d'%e if e>1 else '') for k,e in m) or '1'
            print('             %-24s %s'%(nm,c if c<10**12 else 'C[%s..]'%str(c)[:10]))
# which unknowns occur ONLY in eq8680 among the linear equations
print('\nunknowns whose linear footprint contains 8680:')
for c,eqs in col.items():
    if 8680 in eqs and len(eqs)<=12: print('   x%-6d %d eqs %s'%(NB[c],len(eqs),eqs))
