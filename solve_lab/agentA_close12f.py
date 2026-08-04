#!/usr/bin/env python3
"""Core-fix with UNPINNED-knob preference. Precompute free inputs that appear in any
difference-check (pinned); prefer knobs NOT pinned (truly-free slack) to avoid cascades.
Anchors: x_3558=0 via x_24908 knob x_31339 (keep x_16742); x_29322=0 via x_14853."""
import json, re, sys
from collections import deque, Counter
from agentA_harness import (p, order, definer, gates, freeinp, backward_cone,
                            load_solution, forward, eqcode, eqvars, lines, NEQ, NVARS)
sys.setrecursionlimit(1000000)
CORE = {2071,4573,7123,7469,11854,13660,15299,16622,17726,21382,22093,25480,25539,28653,29437,31061,32894,32916,34517,34892}
gvids = {t: gates[definer[t]][2] for t in order}
base = load_solution('best/new_instance_partial_39013.json'); forward(base)
def count(v):
    ns={'__builtins__':{},'v':v}; return set(i for i in range(NEQ) if eval(eqcode[i],ns)!=0)
allvarcount = Counter()
for i in range(NEQ):
    for x in eqvars[i]: allvarcount[x]+=1
# pinned free inputs = appear in ANY difference (x_a)-(x_b), and const-pinned loads
const_pinned=set(); pinned_free=set()
_dp=re.compile(r'\(x_(\d+)\)\s*-\s*\(x_(\d+)\)')
_bp=re.compile(r'\(x_(\d+)\)\s*-\s*\(?-?\d{20,}'); _bp2=re.compile(r'-?\d{20,}\)?\s*-\s*\(x_(\d+)\)')
for i in range(NEQ):
    L=lines[i]
    for m in _dp.finditer(L):
        a,b=int(m.group(1)),int(m.group(2))
        if a in freeinp: pinned_free.add(a)
        if b in freeinp: pinned_free.add(b)
    if len(L)<4000:
        for m in _bp.finditer(L): const_pinned.add(int(m.group(1)))
        for m in _bp2.finditer(L): const_pinned.add(int(m.group(1)))
print(f"pinned_free={len(pinned_free&freeinp)} const_pinned={len(const_pinned&freeinp)}")
_,frS=backward_cone(35389); _,frT=backward_cone(6671)
QUOT={30317,2936,5146}
protect=(set(frS)|set(frT)|QUOT|const_pinned|{22152,33462,16742,12186,24908,3558,29322})
def free_cone(r):
    seen=set(); st=[r]
    while st:
        u=st.pop()
        if u in seen: continue
        seen.add(u)
        for w in gvids.get(u,()):
            if w not in seen: st.append(w)
    return [u for u in seen if u in freeinp]
_kc={}
def additive_knob(v,gate):
    if gate in _kc: return _kc[gate]
    g0=v[gate]; cands=[]
    for w in free_cone(gate):
        if w in protect: continue
        old=v[w]; v[w]=old+1; forward(v); dd=v[gate]-g0; v[w]=old
        if dd in (1,-1):
            cands.append((1 if w in pinned_free else 0, allvarcount[w], w, dd))  # unpinned first
    forward(v); cands.sort()
    res=(cands[0][2],cands[0][3]) if cands else (None,None)
    _kc[gate]=res; return res
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

v=base[:]
v[31339]=base[31339]+(base[16742]-base[24908]); forward(v)  # x_24908 -> x_16742, x_3558=0
v[14853]=v[12186]; forward(v); set_quot(v)                  # x_29322=0
protect|={31339,14853}
cur=count(v); print(f"after anchors: {NEQ-len(cur)}/{NEQ} ({len(cur)} fail) core={len(cur&CORE)} S,T0={v[35389]%p==0},{v[6671]%p==0}")
dirty=[i for i in range(NVARS) if i in freeinp and v[i]!=base[i] and i not in protect]
queue=deque(dirty); assigned=set(dirty); aligned=set(); best=(len(cur),[x for x in v]); rounds=0
while queue and rounds<40000:
    rounds+=1
    f=queue.popleft(); fval=v[f]
    partners=set()
    for i in eqs_with(f): partners|=partners_in(i,f)
    for G in partners:
        if G in freeinp or G in aligned: continue
        if v[G]==fval: continue
        w,dd=additive_knob(v,G)
        if w is None: continue
        old=v[w]; v[w]=old+dd*(fval-v[G]); forward(v); set_quot(v)
        if v[G]!=fval: v[w]=old; forward(v); set_quot(v); continue
        aligned.add(G); cur=count(v)
        if w not in assigned: assigned.add(w); queue.append(w)
        if len(cur)<best[0]: best=(len(cur),[x for x in v])
    if rounds%50==0 or len(cur)<8:
        print(f" round {rounds}: {NEQ-len(cur)}/{NEQ} ({len(cur)} fail) best={NEQ-best[0]} q={len(queue)} assigned={len(assigned)}",flush=True)
    if not cur: break
v=best[1]; cur=count(v)
print(f"FINAL: {NEQ-len(cur)}/{NEQ} ({len(cur)} fail) core={len(cur&CORE)} S,T0={v[35389]%p==0},{v[6671]%p==0}: {sorted(cur)}")
if NEQ-len(cur)>39022:
    json.dump({f"x_{i}":v[i] for i in range(NVARS) if v[i]!=0}, open(f'best_agentA_{NEQ-len(cur)}.json','w'))
    print(f"SAVED best_agentA_{NEQ-len(cur)}.json")
