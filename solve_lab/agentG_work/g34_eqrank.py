"""Equation-space check: the residual atoms are polynomials in (A,B).  Collect every
equation containing any of them, express its atom combination as alpha*A + beta*B, and
compute the rank of {(alpha,beta)}.  Rank 2 => A=B=0 is FORCED at the equation level."""
import os, sys, pickle, collections
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentG_work')
import g29_frame as F
import gsym2 as G
from gsym2 import L, ad, P
arg=sys.argv[1] if len(sys.argv)>1 else '-'
FL=[int(x) for x in arg.split(',') if x] if arg!='-' else []
v=list(F.v0)
for b in FL: v[b]=1-v[b]
ad.fwd(v,rounds=8)
val,sk=G.build(v,F.NB,cap=6)
# symbolic value of EVERY atom (gates too) so equations can be evaluated symbolically
allsym={}
for a in range(L.NA):
    allsym[a]=G.evalatom(a,val,6)
# residual (nonzero) atoms
nz=[a for a in range(L.NA) if not (isinstance(allsym[a],int) and allsym[a]%P==0)]
print('atoms with nonzero-or-symbolic value: %d'%len(nz))
# equations whose symbolic atom combination is not identically zero
bad=[]
for i,(m,sq,co) in enumerate(L.eq_atoms):
    s=0
    for a,c in co.items():
        t=allsym[a]
        if isinstance(t,int):
            if t%P: s=G.padd(s,(c*t)%P)
        else:
            s=G.padd(s,{mm:(c*cc)%P for mm,cc in t.items()})
    if not (isinstance(s,int) and s%P==0):
        bad.append((i,s))
print('equations whose symbolic value is NOT identically zero: %d'%len(bad))
consts=[(i,s%P) for i,s in bad if isinstance(s,int)]
print('  of which unconditional nonzero constants: %d'%len(consts))
polys=[(i,s) for i,s in bad if not isinstance(s,int)]
print('  of which symbolic: %d'%len(polys))
vars_=sorted({F.NB[k] for i,s in polys for mm in s for k,_ in mm})
print('  variables involved:',vars_[:20],'...' if len(vars_)>20 else '')
pickle.dump({'bad':bad,'NB':F.NB,'flip':FL},open('eqres_%s.pkl'%(arg.replace(',','_')),'wb'))
