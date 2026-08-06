#!/usr/bin/env python3
"""Bounded JOINT Newton: from agentA_39021+loads=CONST, each step target ALL current non-core breaks
jointly (rebuild inner-root Jacobian over closure handles, Dixon pivot-solve, apply, forward). Track
best. Different from greedy one-at-a-time partner alignment. Save if > 39022."""
import json, re, ast, sys, time
from collections import defaultdict
from agentA_harness import (p, order, definer, gates, freeinp, load_solution, forward,
                            eqcode, eqvars, lines, NEQ, NVARS)
sys.setrecursionlimit(1000000)
CONST1=97171863764434070215824145711260403004952728652948669662983319257693684265837195009100680
CONST2=126767545623909574255290391153759363968073470399639361054829680359428658595949132261910506
CORE={2071,4573,7123,7469,11854,13660,15299,16622,17726,21382,22093,25480,25539,28653,29437,31061,32894,32916,34517,34892}
QUOT={30317,2936,5146}; PROT={22152,33462}|QUOT|{12186,14853,16742}
gvids={t:gates[definer[t]][2] for t in order}
ns={'__builtins__':{}}; VAR=re.compile(r'x_(\d+)')
def count(v):
    ns['v']=v; return set(i for i in range(NEQ) if eval(eqcode[i],ns)!=0)
def set_quot(v):
    if v[11150]%p==0: v[30317]=-(v[11150])//p
    if (537773*v[37758])%p==0: v[2936]=(537773*v[37758])//p
    if v[25739]%(6672769*p)==0: v[5146]=v[25739]//(6672769*p)
def inner_src(lhs):
    node=ast.parse(lhs,mode='eval').body
    while isinstance(node,ast.BinOp) and isinstance(node.op,ast.Mult):
        a,b=node.left,node.right
        ca=isinstance(a,ast.Constant) or (isinstance(a,ast.UnaryOp) and isinstance(a.operand,ast.Constant))
        cb=isinstance(b,ast.Constant) or (isinstance(b,ast.UnaryOp) and isinstance(b.operand,ast.Constant))
        if ca and not cb: node=b
        elif cb and not ca: node=a
        elif ast.unparse(a)==ast.unparse(b): node=a
        else: break
    return node
_rc={}
def rootcode(i):
    if i not in _rc:
        _rc[i]=compile(VAR.sub(r'v[\1]',ast.unparse(inner_src(lines[i].rsplit('=',1)[0]))),'<e>','eval')
    return _rc[i]
gcode_k=[compile(VAR.sub(r'v[\1]',gates[definer[order[k]]][1]),'<r>','eval') for k in range(len(order))]
consumers=defaultdict(list)
for k,t in enumerate(order):
    for u in gvids[t]: consumers[u].append(k)
_dk={}
def downstream_ks(w):
    if w in _dk: return _dk[w]
    aff=set(); st=[w]; seen=set()
    while st:
        x=st.pop()
        if x in seen: continue
        seen.add(x)
        for k in consumers.get(x,()):
            if k not in aff: aff.add(k); st.append(order[k])
    _dk[w]=sorted(aff); return _dk[w]
def partial_forward(v,ks):
    ns['v']=v
    for k in ks: v[order[k]]=eval(gcode_k[k],ns)
inv=lambda a: pow(a%p,p-2,p)
eqbyvar=defaultdict(set)
for i in range(NEQ):
    for var in eqvars[i]: eqbyvar[var].add(i)
def matinv(Mmat):
    r=len(Mmat); Aug=[[Mmat[i][j]%p for j in range(r)]+[1 if j==i else 0 for j in range(r)] for i in range(r)]
    for c in range(r):
        piv=None
        for i in range(c,r):
            if Aug[i][c]%p: piv=i;break
        if piv is None: return None
        Aug[c],Aug[piv]=Aug[piv],Aug[c]; ivv=inv(Aug[c][c]); Aug[c]=[(x*ivv)%p for x in Aug[c]]
        for i in range(r):
            if i!=c and Aug[i][c]%p:
                f=Aug[i][c]; Aug[i]=[(Aug[i][k]-f*Aug[c][k])%p for k in range(2*r)]
    return [[Aug[i][r+j] for j in range(r)] for i in range(r)]
def dixon(Mmat,b,steps=16):
    r=len(Mmat); Mi=matinv(Mmat)
    if Mi is None: return None
    x=[0]*r; bb=b[:]; mod=1
    for _ in range(steps):
        bm=[bb[i]%p for i in range(r)]
        xi=[sum(Mi[i][k]*bm[k] for k in range(r))%p for i in range(r)]
        for i in range(r): x[i]+=mod*xi[i]
        nb=[]
        for i in range(r):
            s=bb[i]-sum(Mmat[i][k]*xi[k] for k in range(r))
            if s%p!=0: return None
            nb.append(s//p)
        bb=nb; mod*=p
        if all(z==0 for z in bb): break
    half=mod//2; y=[]
    for xi in x:
        xi%=mod
        if xi>half: xi-=mod
        y.append(xi)
    return y
def gfp_pivots(Jm,H):
    m=len(Jm); n=len(H)
    Mx=[row[:] for row in Jm]; rowmap=list(range(m)); pr=0; pivrows=[]; pivcols=[]
    for c in range(n):
        piv=None
        for i in range(pr,m):
            if Mx[i][c]%p: piv=i;break
        if piv is None: continue
        Mx[pr],Mx[piv]=Mx[piv],Mx[pr]; rowmap[pr],rowmap[piv]=rowmap[piv],rowmap[pr]
        ivv=inv(Mx[pr][c]); Mx[pr]=[(x*ivv)%p for x in Mx[pr]]
        for i in range(m):
            if i!=pr and Mx[i][c]%p:
                f=Mx[i][c]; Mx[i]=[(Mx[i][k]-f*Mx[pr][k])%p for k in range(n)]
        pivrows.append(rowmap[pr]); pivcols.append(c); pr+=1
        if pr>=m: break
    return pivrows,pivcols

v=load_solution('best_agentA_39021.json'); v[22152]=CONST2; v[33462]=CONST1; forward(v); set_quot(v)
best_fail=None; best_v=None
for it in range(10):
    ns['v']=v; F=count(v); nc=sorted(F-CORE); core=F&CORE
    if best_fail is None or len(F)<best_fail:
        best_fail=len(F); best_v=v[:]
    print(f"iter {it}: {NEQ-len(F)}/{NEQ} ({len(F)} fail) core={len(core)} noncore={len(nc)}", flush=True)
    if not nc and not core: break
    target=nc if nc else sorted(core)
    # handles: free inputs in target eqs (1-hop closure), minus PROT
    H=set()
    for i in target: H|=(eqvars[i]&freeinp)
    cons=set(target)
    for h in list(H): cons|=eqbyvar.get(h,set())
    for i in cons: H|=(eqvars[i]&freeinp)
    H=sorted(H-PROT)
    base={i:eval(rootcode(i),ns)%p for i in target}
    dks={h:downstream_ks(h) for h in H}
    Jrows=[[0]*len(H) for _ in target]; tidx={t:i for i,t in enumerate(target)}
    tset=set(target)
    for hj,h in enumerate(H):
        o=v[h]; v[h]=o+1; partial_forward(v,dks[h]); ns['v']=v
        aff=set(eqbyvar.get(h,()))&tset
        for k in dks[h]: aff|=(eqbyvar.get(order[k],set())&tset)
        for i in aff:
            dd=(eval(rootcode(i),ns)-base[i])%p
            if dd: Jrows[tidx[i]][hj]=dd
        v[h]=o; partial_forward(v,dks[h]); ns['v']=v
    pr,pc=gfp_pivots(Jrows,H); rank=len(pr)
    Msq=[[Jrows[pr[i]][pc[j]] for j in range(rank)] for i in range(rank)]
    rhs=[(-base[target[pr[i]]])%p for i in range(rank)]
    y=dixon(Msq,rhs)
    if y is None: print("  dixon failed"); break
    for j in range(rank): v[H[pc[j]]]+=y[j]
    forward(v); set_quot(v)
v=best_v; ns['v']=v; F=count(v)
print(f"BEST: {NEQ-len(F)}/{NEQ} ({len(F)} fail) core={len(F&CORE)} noncore={sorted(F-CORE)}")
if len(F)==0:
    json.dump({f"x_{i}":v[i] for i in range(NVARS) if v[i]!=0}, open('best_agentC_39033.json','w')); print("*** FULL WIN ***")
elif NEQ-len(F)>39022:
    json.dump({f"x_{i}":v[i] for i in range(NVARS) if v[i]!=0}, open(f'best_agentC_{NEQ-len(F)}.json','w')); print(f"IMPROVED best_agentC_{NEQ-len(F)}.json")
