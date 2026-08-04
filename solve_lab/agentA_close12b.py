#!/usr/bin/env python3
"""From 39021 + loads set to CONST1/CONST2 (fixes 12, breaks 24), re-run the partner-
alignment cascade using ONLY clean (non-const-pinned) knobs, to re-align the 24."""
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

def count(v):
    ns={'__builtins__':{},'v':v}; return set(i for i in range(NEQ) if eval(eqcode[i],ns)!=0)

# protect: S,T criticals + the two loads (now at CONST) + const-pinned message endpoints
crit = {16742,14853,12186,24908,22152,33462}
for r in [35389,6671]:
    _,fr=backward_cone(r); crit|=set(fr)
QUOT={30317,2936,5146}
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
protect=set(crit)|QUOT|const_pinned

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
def additive_knob(v,gate):
    if gate in _knob: return _knob[gate]
    g0=v[gate]; cands=[]
    for w in free_cone(gate):
        if w in protect: continue
        old=v[w]; v[w]=old+1; forward(v); dd=v[gate]-g0; v[w]=old
        if dd in (1,-1): cands.append((allvarcount[w],w,dd))
    forward(v); cands.sort()
    res=(cands[0][1],cands[0][2]) if cands else (None,None)
    _knob[gate]=res; return res
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

# start: 39021 + loads to CONST
v=load_solution('best_agentA_39021.json'); v[22152]=CONST2; v[33462]=CONST1
forward(v); set_quot(v)
cur=count(v)
print(f"start (39021+loads=CONST): {NEQ-len(cur)}/{NEQ} ({len(cur)} fail) core={len(cur&CORE)}")

# dirty free inputs = differ from baseline; seed queue
dirty=[i for i in range(NVARS) if i in freeinp and v[i]!=base[i]]
queue=deque(dirty)
assigned={i:v[i] for i in dirty}
aligned_gates=set(); best=(len(cur),[x for x in v]); rounds=0
while queue and rounds<40000:
    rounds+=1
    f=queue.popleft(); fval=v[f]
    partners=set()
    for i in eqs_with(f): partners|=partners_in(i,f)
    for G in partners:
        if G in freeinp or G in aligned_gates: continue
        if v[G]==fval: continue
        w,dd=additive_knob(v,G)
        if w is None: continue
        old=v[w]; v[w]=old+dd*(fval-v[G]); forward(v); set_quot(v)
        if v[G]!=fval:
            v[w]=old; forward(v); set_quot(v); continue
        aligned_gates.add(G); cur=count(v)
        if w not in assigned: assigned[w]=v[w]; queue.append(w)
        if len(cur)<best[0]: best=(len(cur),[x for x in v])
    if rounds%50==0 or len(cur)<10:
        print(f" round {rounds}: {NEQ-len(cur)}/{NEQ} ({len(cur)} fail) best={NEQ-best[0]} q={len(queue)} assigned={len(assigned)}",flush=True)
    if not cur:
        print("ALL SOLVED!"); best=(0,[x for x in v]); break
v=best[1]; cur=count(v)
print(f"FINAL: {NEQ-len(cur)}/{NEQ} ({len(cur)} fail) core={len(cur&CORE)}: {sorted(cur)}")
if NEQ-len(cur)>39021:
    json.dump({f"x_{i}":v[i] for i in range(NVARS) if v[i]!=0}, open(f'best_agentA_{NEQ-len(cur)}.json','w'))
    print(f"SAVED best_agentA_{NEQ-len(cur)}.json")
