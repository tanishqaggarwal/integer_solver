"""Widen the departure support: report per-unknown linear cost for the unknowns that
appear in the region's equations but were excluded by the 'footprint inside' test."""
import os, sys, pickle, collections
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentG_work')
import gsym2 as G
from gsym2 import L, ad, P
D=pickle.load(open('/home/user/integer_solver/solve_lab/agentG_work/coset_model.pkl','rb'))
Lin=pickle.load(open('/home/user/integer_solver/solve_lab/agentG_work/coset_lin.pkl','rb'))
NB=D['NB']; lin=D['lin']; non=D['non']; x0=Lin['x0']; col=Lin['col']
n=len(NB); ix={u:i for i,u in enumerate(NB)}
# is eq8680 a perfect square?
f=dict(non)[8680] if 8680 in dict(non) else None
mons=sorted(f)
a=f.get(((ix[3629],2),),0); b=f.get(((ix[3629],1),(ix[8976],1)),0); c=f.get(((ix[8976],2),),0)
disc=(b*b-4*a*c)%P
print('eq8680 = %d*x3629^2 + C*x3629*x8976 + C*x8976^2 ; discriminant mod p = %d'%(a,disc))
print('  -> PERFECT SQUARE (double root)' if disc==0 else '  -> not a perfect square')
if disc==0:
    r=(-b)*pow(2*a,-1,P)%P
    print('  eq8680 = %d*(x3629 - %d*x8976)^2 ; so eq8680=0 pins x3629 = %d * x8976'%(a,r,r))
S17=[int(x) for x in open('/home/user/integer_solver/solve_lab/agentG_work/extsup.txt').read().split(',')]
print('\nunknowns appearing in the twelve region equations but NOT in the closed support:')
REG=[2554,6816,8124,9123,9421,12231,12270,12350,14584,18673,22044,29125,8680]
inreg=set()
allq={i:f for i,f in lin}; allq.update({i:f for i,f in non})
for q in REG:
    f=allq.get(q)
    if f is None or isinstance(f,int): continue
    inreg |= {NB[k] for m in f for k,_ in m}
for u in sorted(inreg):
    tag='' if u in S17 else '   <-- OUTSIDE the closed support'
    print('   x%-6d  %d linear equations%s'%(u,len(col.get(ix[u],[])),tag))
new=[u for u in sorted(inreg) if u not in S17]
SUP=sorted(set(S17)|set(new))
print('\nwidened support (%d unknowns): %s'%(len(SUP),','.join(map(str,SUP))))
open('/home/user/integer_solver/solve_lab/agentG_work/sup_tierA.txt','w').write(','.join(map(str,SUP)))
for u in new: print('   new x%-6d footprint %s'%(u,col.get(ix[u],[])[:20]))
