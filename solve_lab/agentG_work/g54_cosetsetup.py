"""Equation-level exact model in the deliverable's frame, dumped for coset decoding.

Frame: base state s10/AG_39013.json (boolean free inputs fixed), with the five gate
outputs 7068, 28730, 29854, 31864, 642 DETACHED (turned into free unknowns).  In that
frame the 39,026 deliverable is a point of the parameter space (verified: 0 atom
mismatches over all 42,267 atoms, exactly 7 nonzero equations).

Produces: every equation's inner sum as an exact polynomial over F_p in the unknowns,
split into linear and nonlinear, plus the deliverable's own coordinate vector.
"""
import os, sys, pickle, time, collections
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentG_work')
import gframe, gGclose, gsolve
import gsym2 as G
from gsym2 import L, ad, P
OUT='/home/user/integer_solver/solve_lab/agentG_work/coset_model.pkl'
DET=[7068,28730,29854,31864,642]
BASE='/home/user/integer_solver/solve_lab/s10/AG_39013.json'
DEL='/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json'
v=L.load(BASE); ad.fwd(v,rounds=6)
vd=L.load(DEL)
FREE=[u for u in range(L.NVARS) if u not in L.definer or u in set(DET)]
NB=sorted({u for u in FREE if not gGclose.isbool(u)} | set(DET))
n=len(NB); ix={u:i for i,u in enumerate(NB)}
print('unknowns (non-boolean free inputs + detached): %d'%n,flush=True)
t0=time.time()
val=gframe.build_detached(v,set(DET),NB,cap=6)
allsym=[G.evalatom(a,val,6) for a in range(L.NA)]
lin=[];non=[]
for i,(m,sq,co) in enumerate(L.eq_atoms):
    s=0
    for a,c in co.items():
        t=allsym[a]
        if isinstance(t,int):
            if t%P: s=G.padd(s,(c*t)%P)
        else: s=G.padd(s,{mm:(c*cc)%P for mm,cc in t.items()})
    if isinstance(s,int):
        if s%P: non.append((i,s%P))
        continue
    (lin if G.deg(s)==1 else non).append((i,s))
print('non-trivial equations %d = linear %d + nonlinear/const %d  (%.0fs)'%(len(lin)+len(non),len(lin),len(non),time.time()-t0),flush=True)
# the deliverable's own point
pt=[vd[u]%P for u in NB]
def ev(f):
    if isinstance(f,int): return f%P
    s=0
    for m,c in f.items():
        t=c
        for k,e in m: t=t*pow(pt[k],e,P)%P
        s=(s+t)%P
    return s
badlin=[i for i,f in lin if ev(f)]
badnon=[i for i,f in non if ev(f)]
print('AT THE DELIVERABLE: violated linear = %s'%badlin)
print('                    violated nonlinear = %s'%badnon)
print('                    total = %d'%(len(badlin)+len(badnon)))
pickle.dump({'NB':NB,'lin':lin,'non':non,'pt':pt,'badlin':badlin,'badnon':badnon,'DET':DET},open(OUT,'wb'))
print('dumped',OUT)
