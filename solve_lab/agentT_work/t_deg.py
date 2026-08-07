#!/usr/bin/env python3
"""AUDIT T24c -- is L's 'degree never exceeds 3' real, or an artifact of fitting with deg=4?
fit() samples t=0..4 -- exactly 5 points -- so it can only SEE degree <= 4.  A true degree-5+
polynomial would be silently aliased to a wrong degree-4 fit.  Re-fit with more points."""
import sys, os, json, time
from math import gcd
L='/home/user/integer_solver/solve_lab/agentL_work'
os.chdir(L); sys.path.insert(0,L)
src=open(os.path.join(L,'mkassign2.py')).read().split('#MAINSTART')[0]
g={'__name__':'mkdriver'}
exec(compile(src,'pre','exec'),g)
p=g['p']; E=g['E']; NV=g['NV']; SL=g['SL']; SHIFT=g['SHIFT']
ORIENT=g['ORIENT']; T1=g['T1']; T2=g['T2']; assignment=g['assignment']; relift=g['relift']
atomvalvars=g['atomvalvars']; vars_of=g['vars_of']
lines=open(os.path.join(L,'solve927.py')).read().split(chr(10))
keep=[]
for ln in lines:
    if ln.startswith('src=open(') or ln.startswith('exec(src)'): continue
    if ln.startswith("if __name__"): break
    keep.append(ln)
exec(compile(chr(10).join(keep),'s927','exec'),g)
probe=g['probe']
def fitn(vv,i,w,deg):
    ys=[]
    for t in range(deg+1):
        y=probe(vv,i,[w],[t])
        if y%p: return None,None
        ys.append(y//p)
    dd=[ys[:]]
    for k in range(deg):
        dd.append([dd[k][j+1]-dd[k][j] for j in range(len(dd[k])-1)])
    co=[dd[k][0] for k in range(deg+1)]
    return co, max([k for k in range(deg+1) if co[k]!=0], default=0)
S=[24601,2081]
v,isl,valn=assignment(set(S),ORIENT); v[24468]=T1; v[18956]=T2
vv=[0]*NV
for k,x in v.items(): vv[k]=x
for rd in range(60):
    bad=relift(vv)
    if not bad: break
    r=E.run(vv); fx=0
    for a in bad:
        i=E.residx[a]; cur=r[i]; sm=abs(SL[a])
        if cur%p: continue
        imm=[q for q in vars_of(E.atoms[a]) if q in SHIFT]
        for w in imm+[q for q in atomvalvars[a] if q in SHIFT and q not in imm]:
            old=vv[w]; vv[w]=old+p; dl=E.run(vv)[i]-cur; vv[w]=old
            if dl==0: continue
            gg=gcd(dl,sm)
            if cur%gg: continue
            mm=sm//gg
            t=(-(cur//gg))*pow((dl//gg)%mm,-1,mm)%mm if mm>1 else 0
            vv[w]=old+p*t; fx+=1; break
    if fx==0: break
left=relift(vv); r=E.run(vv)
stuck=[a for a in left if r[E.residx[a]]%p==0]
print('stuck conditions at the greedy fixpoint: %d'%len(stuck),flush=True)
for a in stuck:
    i=E.residx[a]; c=abs(SL[a])//p
    ws=[q for q in vars_of(E.atoms[a]) if q in SHIFT]
    ws+=[q for q in atomvalvars[a] if q in SHIFT and q not in ws]
    print('\n  condition c=%d  atom %s'%(c,a[:60]),flush=True)
    for w in ws[:6]:
        row=[]
        for D in (4,6,8,10):
            co,td=fitn(vv,i,w,D)
            row.append('deg<=%d -> top %s'%(D,td if co is not None else 'n/a'))
        print('     wire x%-6d  %s'%(w,' | '.join(row)),flush=True)
