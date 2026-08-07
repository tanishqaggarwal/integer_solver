"""Exact max-satisfy over the enlarged residual region.
Rows: equations of R*, affine in the knob vector delta (absolute knob values).
Solve: which subsets of rows are simultaneously satisfiable over Z?"""
import sys, collections, itertools, json; sys.path.insert(0,'.')
import env, lib as L
from amk_model import build, knobpoly, v0, av0, E
from fractions import Fraction
P=env.P

EXTRA=[37887,41906]
A,K,R=build(EXTRA)
Aset=set(A)
print('atoms=%d knobs=%d eqs=%d'%(len(A),len(K),len(R)))
kp={a:knobpoly(a,K,v0) for a in A}
QUAD=[a for a in A if any(len(m)>1 for m in kp[a])]
print('atoms nonlinear in knobs:',QUAD)

# affine part per atom: const + sum c_i * d_i     (only for non-quad atoms)
def affine(a):
    P0=kp[a]; c=P0.get((),0); lin=[0]*len(K)
    for m,cc in P0.items():
        if len(m)==1: lin[m[0]]=cc
    return c,lin

rows=[]   # (eq, const, lin[], has_quad_atom)
for e in R:
    m,sq,co=L.eq_atoms[e]
    # foreign atoms (not in A) must be zero and stay zero
    for a in co:
        if a not in Aset:
            assert av0[a]==0, (e,a)
    c=0; lin=[0]*len(K); hq=False
    for a,cc in co.items():
        if a not in Aset: continue
        if a in QUAD: hq=True; continue
        c0,l0=affine(a); c+=cc*c0
        for i in range(len(K)): lin[i]+=cc*l0[i]
    rows.append({'eq':e,'c':c,'lin':lin,'quad':hq})
d0=[v0[u] for u in K]
for r in rows:
    r['cur']=r['c']+sum(a*b for a,b in zip(r['lin'],d0))
bad=[r['eq'] for r in rows if r['cur']!=0 and not r['quad']]
print('rows=%d  currently violated (affine part):'%len(rows), bad)
print('rows with quad atom:',[r['eq'] for r in rows if r['quad']])
json.dump({'K':K,'A':A,'R':R,
           'rows':[{'eq':r['eq'],'c':str(r['c']),'lin':[str(x) for x in r['lin']],'quad':r['quad']} for r in rows],
           'd0':[str(x) for x in d0]}, open('/home/user/integer_solver/solve_lab/agentA_work/model22.json','w'))
print('saved model22.json')
