"""Is the 39,026 deliverable inside my symbolic parameter space?"""
import os, sys
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentG_work')
import gframe, gGclose
import gsym2 as G
from gsym2 import L, ad, P
DEL='/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json'
BASE='/home/user/integer_solver/solve_lab/s10/AG_39013.json'
vd=L.load(DEL); vb=L.load(BASE); ad.fwd(vb,rounds=6)
avd=L.all_atom_values(vd)
nz=[a for a in range(L.NA) if avd[a]]
print('deliverable nonzero atoms:',nz,'failing',len(L.failing_eqs(avd)))
FREE=[u for u in range(L.NVARS) if u not in L.definer]
BOOL=[u for u in FREE if gGclose.isbool(u)]
diffb=[u for u in BOOL if (vd[u]%P)!=(vb[u]%P)]
print('boolean free inputs differing deliverable vs AG_39013:',len(diffb),diffb[:20])
print('booleans set to 1 in deliverable:',[u for u in BOOL if vd[u]==1][:20])
print('booleans set to 1 in AG_39013   :',[u for u in BOOL if vb[u]==1][:20])
DET=[L.atom_out[a][1] for a in nz if a in L.atom_out]
print('detach set from deliverable:',DET)
# evaluate my symbolic system at the deliverable's point
NB=[u for u in FREE if not gGclose.isbool(u)]+[u for u in DET if u not in FREE]
NB=sorted(set(NB))
val=gframe.build_detached(vd, set(DET), NB, cap=6)
allsym=[G.evalatom(a,val,6) for a in range(L.NA)]
pt={i:(vd[u]%P) for i,u in enumerate(NB)}
def ev(f):
    if isinstance(f,int): return f%P
    s=0
    for m,c in f.items():
        t=c
        for k,e in m: t=t*pow(pt[k],e,P)%P
        s=(s+t)%P
    return s
bad=0; mism=0
for a in range(L.NA):
    if ev(allsym[a])!=avd[a]%P:
        mism+=1
        if mism<5: print('  atom mismatch a%d'%a)
print('atom-value mismatches between symbolic model at deliverable point and truth:',mism)
cnt=0
for i,(m,sq,co) in enumerate(L.eq_atoms):
    s=0
    for a,c in co.items():
        t=allsym[a]
        if isinstance(t,int): s=(s+c*t)%P
        else: s=(s+c*ev(t))%P
    if s%P: cnt+=1
print('equations nonzero mod p at the deliverable point (my model): %d'%cnt)
