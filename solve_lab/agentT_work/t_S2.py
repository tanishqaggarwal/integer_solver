#!/usr/bin/env python3
"""AUDIT T24 -- agent L's |S|=2 closure over Z.

L reports: for ON = {24601, 2081}, all 927 c>1 conditions discharged, '0 undischarged, 2 nonzero
atoms = the target congruences'.  The whole integer half of the reduction rests on it.

THE PREMISE NOBODY TESTED: every number in that result is computed inside L's OWN engine
(`E.run`, a 9,032-residual-atom model).  `solve927.py` prints and **dumps no assignment**, so the
closure has never been put in front of `checker.py`.  L's |S|=1 result WAS checker-verified
(assign_L1.json, 15 failing); the |S|=2 closure was not, and assign_L2.json predates it by 2h.

Reproduce the closure, dump the assignment, and check it against the real instance.
Read-only w.r.t. agentL_work (PYTHONDONTWRITEBYTECODE=1; all output into agentT_work)."""
import sys, os, json, time
from math import gcd
L='/home/user/integer_solver/solve_lab/agentL_work'
T='/home/user/integer_solver/solve_lab/agentT_work'
os.chdir(L); sys.path.insert(0,L)
src=open(os.path.join(L,'mkassign2.py')).read().split('#MAINSTART')[0]
g={'__name__':'mkdriver'}
exec(compile(src,'mkassign2_prefix','exec'),g)
p=g['p']; E=g['E']; NV=g['NV']; M=g['M']; SL=g['SL']; SHIFT=g['SHIFT']
ORIENT=g['ORIENT']; T1=g['T1']; T2=g['T2']; assignment=g['assignment']; relift=g['relift']
atomvalvars=g['atomvalvars']; vars_of=g['vars_of']
lines=open(os.path.join(L,'solve927.py')).read().split(chr(10))
keep=[]
for ln in lines:
    if ln.startswith('src=open(') or ln.startswith('exec(src)'): continue
    if ln.startswith("if __name__"): break
    keep.append(ln)
exec(compile(chr(10).join(keep),'solve927_prefix','exec'),g)
factor=g['factor']; probe=g['probe']; fit=g['fit']; solve_one=g['solve_one']
S=[24601,2081]
v,isl,valn=assignment(set(S),ORIENT); v[24468]=T1; v[18956]=T2
vv=[0]*NV
for k,x in v.items(): vv[k]=x
t0=time.time()
for rd in range(60):
    bad=relift(vv)
    if not bad: break
    r=E.run(vv); fx=0
    for a in bad:
        i=E.residx[a]; cur=r[i]; sm=abs(SL[a])
        if cur%p: continue
        imm=[q for q in vars_of(E.atoms[a]) if q in SHIFT]
        for w in imm+[q for q in atomvalvars[a] if q in SHIFT and q not in imm]:
            old=vv[w]; vv[w]=old+p; d=E.run(vv)[i]-cur; vv[w]=old
            if d==0: continue
            gg=gcd(d,sm)
            if cur%gg: continue
            mm=sm//gg
            t=(-(cur//gg))*pow((d//gg)%mm,-1,mm)%mm if mm>1 else 0
            vv[w]=old+p*t; fx+=1; break
    if fx==0: break
print('greedy fixpoint reached (%.0fs)'%(time.time()-t0),flush=True)
for outer in range(12):
    left=relift(vv); r=E.run(vv)
    stuck=[a for a in left if r[E.residx[a]]%p==0]
    if not stuck:
        print('ALL c>1 CONDITIONS DISCHARGED after %d outer rounds'%outer,flush=True); break
    print('round %d: %d stuck'%(outer,len(stuck)),flush=True)
    prog=0
    for a in stuck:
        res=solve_one(vv,a)
        if res and len(res)==3:
            w,t,td=res; vv[w]+=p*t; prog+=1
            print('   SOLVED c=%-9d deg=%d wire x%-6d t=%-10d'%(abs(SL[a])//p,td,w,t),flush=True)
    if prog==0: break
left=relift(vv); r=E.run(vv)
stuck=[a for a in left if r[E.residx[a]]%p==0]
nz=[a for a in E.residx if r[E.residx[a]]]
print('\nIN L\'s MODEL: %d undischarged, %d nonzero atoms of %d'%(len(stuck),len(nz),len(E.residx)),flush=True)
for a in nz: print('    %s'%a[:88],flush=True)
# --- AUDIT ADDITION 1: how many atoms in `left` were EXCLUDED from the stuck count? ---
notmodp=[a for a in left if r[E.residx[a]]%p!=0]
print('\nbad-list entries EXCLUDED from the "stuck" count (residual NOT 0 mod p): %d'%len(notmodp),flush=True)
# --- AUDIT ADDITION 2: dump and check against the REAL instance ---
out=os.path.join(T,'t_S2_assign.json')
json.dump({'x_%d'%i:str(vv[i]) for i in range(NV) if vv[i]}, open(out,'w'))
print('\nwrote %s'%out,flush=True)
