"""Measure, by construction, the true dimension of the reachable equation-residual space.

For every perturbation that leaves the atom support S unchanged (the 'harmless' set from
the census, plus every variable of S and its knobs), record the induced change of the 12
equations of E(S) as a vector in Q^12.  The rank of the collected vectors is the number
of equations we can steer => failing >= |E| - rank.  This is the balance law's dof,
measured rather than assumed.
"""
import os, sys, collections, json, time, random
from fractions import Fraction
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L
BROKEN_GATES={22229,22230,35758,35761,35762}
BASE={22229,22230,35758,35759,35760,35761,35762}
E0=sorted(set().union(*[set(L.atom2eq[a]) for a in BASE]))
v0=L.load(os.path.join(LAB,'best','new_instance_partial_39026.json'))
base_av=L.all_atom_values(v0)
def eqvec(av): return [L.eq_value(e,av) for e in E0]
r0=eqvec(base_av)
print('E0',E0)
print('r0 nonzero at',[E0[i] for i,x in enumerate(r0) if x])

cen={}
for f in ('pa_census_a.json','pa_census_b.json'):
    cen.update(json.load(open(os.path.join(HERE,f))))
harmless=[int(t) for t,r in cen.items() if r['fail']==7 and set(r['nz'])==BASE]
VS=sorted(set(harmless)|{u for a in BASE for u in L.avars[a]})
print('candidate knob variables:',len(VS))

rows=[]; owners=[]
random.seed(7)
DS=[1, 3, random.randrange(1,2**64), random.randrange(1,2**200)]
t0=time.time()
for i,t in enumerate(VS):
    for d in DS:
        w=list(v0)
        try:
            ch,_=L.ripple(w,{t:v0[t]+d},maxsteps=60000,block=BROKEN_GATES)
        except Exception: continue
        cand=set()
        for u in ch: cand.update(L.var_atoms[u])
        nz=set(BASE); ok=True
        av=dict()
        for a in cand:
            nv=L.evalpoly(L.polys[a],w)
            av[a]=nv
            if nv and a not in BASE: ok=False; break
            if not nv: nz.discard(a)
        if not ok: continue
        full=list(base_av)
        for a,x in av.items(): full[a]=x
        vec=[a-b for a,b in zip(eqvec(full),r0)]
        if any(vec):
            rows.append(vec); owners.append((t,d))
    if i%300==0: print(' ',i,len(rows),f'{time.time()-t0:.0f}s',flush=True)
print('collected',len(rows),'nonzero direction vectors')

def rank_and_solve(rows,target):
    n=len(E0)
    M=[[Fraction(x) for x in r] for r in rows]
    piv=[]; used=[]
    R=[]
    for r in M:
        cur=list(r)
        for (pc,pr) in R:
            if cur[pc]:
                f=cur[pc]/pr[pc]
                cur=[a-f*b for a,b in zip(cur,pr)]
        nz=[j for j in range(n) if cur[j]]
        if nz:
            R.append((nz[0],cur))
    rk=len(R)
    # is -r0 in the span?
    cur=[Fraction(-x) for x in target]
    for (pc,pr) in R:
        if cur[pc]:
            f=cur[pc]/pr[pc]
            cur=[a-f*b for a,b in zip(cur,pr)]
    resid=[j for j in range(n) if cur[j]]
    return rk,resid,R

rk,resid,R=rank_and_solve(rows,r0)
print(f'RANK of reachable equation-direction space = {rk}  (of |E|={len(E0)})')
print(f'  => min failing (linear/generic) = |E| - rank = {len(E0)-rk}')
print(f'  residual after projection: nonzero coords {resid} -> r0 in span: {not resid}')
# which owners contribute independent directions
seen=[]; contrib=[]
R2=[]
for vec,ow in zip(rows,owners):
    cur=[Fraction(x) for x in vec]
    for (pc,pr) in R2:
        if cur[pc]:
            f=cur[pc]/pr[pc]; cur=[a-f*b for a,b in zip(cur,pr)]
    nz=[j for j in range(len(E0)) if cur[j]]
    if nz:
        R2.append((nz[0],cur)); contrib.append(ow)
print('independent knobs:',contrib)
json.dump({'rank':rk,'contrib':[[t,str(d)] for t,d in contrib],'E0':E0},open(os.path.join(HERE,'pa_dof.json'),'w'))
print(f'{time.time()-t0:.0f}s done')
