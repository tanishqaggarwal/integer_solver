#!/usr/bin/env python3
"""Clean redo: from best_39013 (loads satisfied, x_33462=CONST1, x_22152=CONST2), apply the regime-1
core move (x_29322->0, x_3558->0), forward, and Dixon-heal the resulting breaks while PROTECTING the
load-critical inputs x_33462,x_22152 (+load cone) and core cone. If consistent -> full solve."""
import json, re, ast, sys, time
from collections import defaultdict
import agentC_common as AC
from agentC_common import (p, order, definer, gates, gcode, forward, val, freeinp, ns, lines,
                           eqcode, eqvars, load_best, CORE, NVARS, rootcode_of, inv, posof,
                           downstream_ks, partial_forward, C1, C2)
CONST1=97171863764434070215824145711260403004952728652948669662983319257693684265837195009100680
CONST2=126767545623909574255290391153759363968073470399639361054829680359428658595949132261910506
QUOT={30317,2936,5146}
gvids={t:gates[definer[t]][2] for t in order}
def free_cone(r):
    seen=set(); st=[r]
    while st:
        u=st.pop()
        if u in seen: continue
        seen.add(u)
        for w in gvids.get(u,()):
            if w not in seen: st.append(w)
    return set(u for u in seen if u in freeinp)
def set_quot():
    if val[11150]%p==0: val[30317]=-(val[11150])//p
    if (537773*val[37758])%p==0: val[2936]=(537773*val[37758])//p
    if val[25739]%(6672769*p)==0: val[5146]=val[25739]//(6672769*p)
def fails():
    ns['v']=val; return set(i for i in range(len(lines)) if eval(eqcode[i],ns)!=0)
eqbyvar=defaultdict(set)
for i in range(len(lines)):
    for var in eqvars[i]: eqbyvar[var].add(i)

load_best(); forward(); ns['v']=val
print(f"best_39013: {len(lines)-len(fails())}/{len(lines)}; x_33462==CONST1:{val[33462]==CONST1}, x_22152==CONST2:{val[22152]==CONST2}")
# regime-1 core move: x_12186 += x_29322 residue -> x_29322=0 ; x_16742 += x_3558 residue -> x_3558=0
da=val[29322]%p; dc=val[3558]%p
val[12186]+=da; val[16742]+=dc; forward(); set_quot(); ns['v']=val
F=fails(); core=F&CORE; nc=F-CORE
print(f"after regime-1 move: {len(lines)-len(F)}/{len(lines)} core={len(core)} noncore={len(nc)}: {sorted(nc)}")
print(f"  S={val[35389]%p}, T={val[6671]%p}")
# PROTECT: loads + their cones + core cone + the moved inputs + quot
loadeqs=[i for i in range(len(lines)) if 33462 in eqvars[i] or 22152 in eqvars[i]]
protect={33462,22152,12186,16742,14853,24908}|QUOT
for i in loadeqs:
    for var in eqvars[i]:
        protect|=free_cone(var)
        if var in freeinp: protect.add(var)
protect|=free_cone(35389)|free_cone(6671)  # core cone
print(f"protect set: {len(protect)} frees (incl load+core cones)")

gcode_k=gcode
def solve_and_apply(target):
    H=set()
    for i in target: H|=(eqvars[i]&freeinp)
    cons=set(target)
    for h in list(H): cons|=eqbyvar.get(h,set())
    for i in cons: H|=(eqvars[i]&freeinp)
    H=sorted(H-protect)
    base={i:eval(rootcode_of(i),ns)%p for i in target}
    dks={h:downstream_ks(h) for h in H}
    Jrows=[[0]*len(H) for _ in target]; tidx={t:i for i,t in enumerate(target)}; tset=set(target)
    for hj,h in enumerate(H):
        o=val[h]; val[h]=o+1; partial_forward(dks[h]); ns['v']=val
        aff=set(eqbyvar.get(h,()))&tset
        for k in dks[h]: aff|=(eqbyvar.get(order[k],set())&tset)
        for i in aff:
            d=(eval(rootcode_of(i),ns)-base[i])%p
            if d: Jrows[tidx[i]][hj]=d
        val[h]=o; partial_forward(dks[h]); ns['v']=val
    # pivots
    m=len(target); n=len(H); Mx=[row[:] for row in Jrows]; rowmap=list(range(m)); pr=0; pivr=[]; pivc=[]
    for c in range(n):
        sel=None
        for i in range(pr,m):
            if Mx[i][c]%p: sel=i;break
        if sel is None: continue
        Mx[pr],Mx[sel]=Mx[sel],Mx[pr]; rowmap[pr],rowmap[sel]=rowmap[sel],rowmap[pr]
        ivv=inv(Mx[pr][c]); Mx[pr]=[(x*ivv)%p for x in Mx[pr]]
        for i in range(m):
            if i!=pr and Mx[i][c]%p:
                f=Mx[i][c]; Mx[i]=[(Mx[i][k]-f*Mx[pr][k])%p for k in range(n)]
        pivr.append(rowmap[pr]); pivc.append(c); pr+=1
        if pr>=m: break
    rank=len(pivr)
    # consistency: all target rows must be in row space (rank == #independent target constraints)
    Msq=[[Jrows[pivr[i]][pivc[j]] for j in range(rank)] for i in range(rank)]
    rhs=[(-base[target[pivr[i]]])%p for i in range(rank)]
    # matinv+dixon
    r=rank; Aug=[[Msq[i][j]%p for j in range(r)]+[1 if j==i else 0 for j in range(r)] for i in range(r)]
    for c in range(r):
        piv=None
        for i in range(c,r):
            if Aug[i][c]%p: piv=i;break
        if piv is None: return None,rank,len(H),False
        Aug[c],Aug[piv]=Aug[piv],Aug[c]; ivv=inv(Aug[c][c]); Aug[c]=[(x*ivv)%p for x in Aug[c]]
        for i in range(r):
            if i!=c and Aug[i][c]%p:
                f=Aug[i][c]; Aug[i]=[(Aug[i][k]-f*Aug[c][k])%p for k in range(2*r)]
    Mi=[[Aug[i][r+j] for j in range(r)] for i in range(r)]
    # check all target rows satisfied by pivot soln (mod p)
    xsol=[sum(Mi[i][k]*rhs[k] for k in range(r))%p for i in range(r)]
    consistent=all(sum(Jrows[ti][pivc[j]]*xsol[j] for j in range(r))%p == (-base[target[ti]])%p for ti in range(m))
    if not consistent: return None,rank,len(H),False
    # dixon lift
    x=[0]*r; bb=rhs[:]; mod=1
    for _ in range(16):
        bm=[bb[i]%p for i in range(r)]
        xi=[sum(Mi[i][k]*bm[k] for k in range(r))%p for i in range(r)]
        for i in range(r): x[i]+=mod*xi[i]
        nb=[]
        for i in range(r):
            s=bb[i]-sum(Msq[i][k]*xi[k] for k in range(r))
            if s%p!=0: return None,rank,len(H),True
            nb.append(s//p)
        bb=nb; mod*=p
        if all(z==0 for z in bb): break
    half=mod//2
    for j in range(r):
        xi=x[j]%mod
        if xi>half: xi-=mod
        val[H[pivc[j]]]+=xi
    forward(); set_quot(); ns['v']=val
    return True,rank,len(H),True

best_state=None; best_fail=len(F)
for it in range(8):
    ns['v']=val; F=fails(); nc=sorted(F-CORE); core=F&CORE
    if len(F)<best_fail: best_fail=len(F); best_state=val[:]
    print(f"heal iter {it}: {len(lines)-len(F)}/{len(lines)} core={len(core)} noncore={len(nc)}", flush=True)
    if not F: break
    target=nc if nc else sorted(core)
    res,rank,nh,cons=solve_and_apply(target)
    print(f"   target={len(target)} handles={nh} rank={rank} consistent={cons} applied={res is not None}", flush=True)
    if res is None:
        # inconsistent - loads protected block the heal; stop
        break
ns['v']=val; F=fails()
if len(F)<best_fail: best_fail=len(F); best_state=val[:]
print(f"FINAL best: {len(lines)-best_fail}/{len(lines)} fail={best_fail}")
if best_state is not None:
    for i in range(NVARS): val[i]=best_state[i]
    ns['v']=val; F=fails()
    if len(F)==0:
        json.dump({f"x_{i}":val[i] for i in range(NVARS) if val[i]!=0}, open('best_agentC_39033.json','w')); print("*** FULL WIN ***")
    elif len(lines)-len(F)>39022:
        json.dump({f"x_{i}":val[i] for i in range(NVARS) if val[i]!=0}, open(f'best_agentC_{len(lines)-len(F)}.json','w')); print(f"IMPROVED -> best_agentC_{len(lines)-len(F)}.json")
    print(f"loads: x_33462==CONST1:{val[33462]==CONST1}, x_22152==CONST2:{val[22152]==CONST2}; core={len(F&CORE)}")
