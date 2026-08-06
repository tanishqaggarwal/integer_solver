#!/usr/bin/env python3
"""Core-fix holding loads fixed. Zero x_3558 via x_24908 knob x_31339 (keeps x_16742 side
=> x_19083, x_33462 untouched). Zero x_29322 via x_12186 or x_14853. Protect ALL const-
pinned load endpoints during cascade heal so x_33462, x_22152 stay at their CONST."""
import json, re, sys
from collections import deque, Counter
from agentA_harness import (p, order, definer, gates, freeinp, backward_cone,
                            load_solution, forward, eqcode, eqvars, lines, NEQ, NVARS)
sys.setrecursionlimit(1000000)
CONST1 = 97171863764434070215824145711260403004952728652948669662983319257693684265837195009100680
CONST2 = 126767545623909574255290391153759363968073470399639361054829680359428658595949132261910506
CORE = {2071,4573,7123,7469,11854,13660,15299,16622,17726,21382,22093,25480,25539,28653,29437,31061,32894,32916,34517,34892}
gvids = {t: gates[definer[t]][2] for t in order}
base = load_solution('best/new_instance_partial_39013.json'); forward(base)
print(f"baseline x_22152==CONST2? {base[22152]==CONST2}  x_33462==CONST1? {base[33462]==CONST1}")
print(f"baseline x_22152={base[22152]}  x_33462={base[33462]}")

def count(v):
    ns={'__builtins__':{},'v':v}; return set(i for i in range(NEQ) if eval(eqcode[i],ns)!=0)

allvarcount = Counter()
for i in range(NEQ):
    for x in eqvars[i]: allvarcount[x]+=1
const_pinned=set()
_bp=re.compile(r'\(x_(\d+)\)\s*-\s*\(?-?\d{20,}'); _bp2=re.compile(r'-?\d{20,}\)?\s*-\s*\(x_(\d+)\)')
for i in range(NEQ):
    L=lines[i]
    if len(L)<4000:
        for m in _bp.finditer(L): const_pinned.add(int(m.group(1)))
        for m in _bp2.finditer(L): const_pinned.add(int(m.group(1)))
print(f"const-pinned loads protected: {len(const_pinned&freeinp)}")

_knob={}
def free_cone(r):
    seen=set(); st=[r]
    while st:
        u=st.pop()
        if u in seen: continue
        seen.add(u)
        for w in gvids.get(u,()):
            if w not in seen: st.append(w)
    return [u for u in seen if u in freeinp]
def additive_knob(v,gate,protect):
    key=gate
    g0=v[gate]; cands=[]
    for w in free_cone(gate):
        if w in protect: continue
        old=v[w]; v[w]=old+1; forward(v); dd=v[gate]-g0; v[w]=old
        if dd in (1,-1): cands.append((allvarcount[w],w,dd))
    forward(v); cands.sort()
    return (cands[0][1],cands[0][2]) if cands else (None,None)
def partners_in(eq_i,f):
    L=lines[eq_i]; out=set()
    for m in re.finditer(r'\(x_'+str(f)+r'\)\s*-\s*\(x_(\d+)\)',L): out.add(int(m.group(1)))
    for m in re.finditer(r'\(x_(\d+)\)\s*-\s*\(x_'+str(f)+r'\)',L): out.add(int(m.group(1)))
    return out
_v2e={}
def eqs_with(var):
    if var not in _v2e: _v2e[var]=[i for i in range(NEQ) if var in eqvars[i]]
    return _v2e[var]
def set_quot(v):
    if v[11150]%p==0: v[30317]=-(v[11150])//p
    if (537773*v[37758])%p==0: v[2936]=(537773*v[37758])//p
    if v[25739]%(6672769*p)==0: v[5146]=v[25739]//(6672769*p)

QUOT={30317,2936,5146}
_,frS=backward_cone(35389); _,frT=backward_cone(6671)
STcone=set(frS)|set(frT)

def run(x3558_side, x29322_side):
    v=base[:]
    # x_3558 = x_24908 - x_16742 -> 0
    if x3558_side=='x24908':   # move x_24908 to x_16742 via x_31339, keep x_16742
        d=1  # x_31339 deriv on x_24908 (assume +1; verify)
        v[31339]=base[31339]+(base[16742]-base[24908]); forward(v)
        anchor3558=[31339]
    else:                       # move x_16742 to x_24908
        v[16742]=base[24908]; forward(v); anchor3558=[16742]
    # x_29322 = x_14853 - x_12186 -> 0
    if x29322_side=='x12186':
        v[12186]=v[14853]; anchor29322=[12186]
    else:
        v[14853]=v[12186]; anchor29322=[14853]
    forward(v); set_quot(v)
    # protect: S,T cone (minus anchors) + QUOT + all const-pinned loads
    protect=(STcone|QUOT|const_pinned) - set(anchor3558) - set(anchor29322)
    protect|= {22152,33462}
    dirty=[i for i in range(NVARS) if i in freeinp and v[i]!=base[i]]
    queue=deque(dirty); assigned={i:v[i] for i in dirty}
    aligned=set(); cur=count(v); best=(len(cur),[x for x in v]); rounds=0
    while queue and rounds<30000:
        rounds+=1
        f=queue.popleft(); fval=v[f]
        partners=set()
        for i in eqs_with(f): partners|=partners_in(i,f)
        for G in partners:
            if G in freeinp or G in aligned: continue
            if v[G]==fval: continue
            w,dd=additive_knob(v,G,protect)
            if w is None: continue
            old=v[w]; v[w]=old+dd*(fval-v[G]); forward(v); set_quot(v)
            if v[G]!=fval: v[w]=old; forward(v); set_quot(v); continue
            aligned.add(G); cur=count(v)
            if w not in assigned: assigned[w]=v[w]; queue.append(w)
            if len(cur)<best[0]: best=(len(cur),[x for x in v])
        if not cur: break
    v=best[1]; cur=count(v)
    stok=(v[35389]%p==0,v[6671]%p==0)
    loadok=(v[22152]==CONST2 and v[33462]==CONST1) if base[22152]==CONST2 else (v[22152]==base[22152] and v[33462]==base[33462])
    return v,cur

for s3,s29 in [('x24908','x12186'),('x24908','x14853')]:
    v,cur=run(s3,s29)
    print(f"\n[{s3} | {s29}]: {NEQ-len(cur)}/{NEQ} ({len(cur)} fail) core={len(cur&CORE)} "
          f"S,T0={v[35389]%p==0},{v[6671]%p==0} x22152ok={v[22152]==base[22152]} x33462ok={v[33462]==base[33462]}")
    print(f"   fails: {sorted(cur)[:20]}")
    if NEQ-len(cur)>39021:
        json.dump({f"x_{i}":v[i] for i in range(NVARS) if v[i]!=0}, open(f'best_agentA_{NEQ-len(cur)}.json','w'))
        print(f"   SAVED best_agentA_{NEQ-len(cur)}.json")
