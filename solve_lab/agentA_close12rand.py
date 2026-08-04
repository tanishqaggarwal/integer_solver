#!/usr/bin/env python3
"""Randomized-knob cascade, multiple seeds. Different knob choices may avoid the
x_4432/x_7068 wall. Anchors: x_3558=0 via x_31339, x_29322=0 via x_14853. Keep best."""
import json, re, sys, random
from collections import deque, Counter
from agentA_harness import (p, order, definer, gates, freeinp, backward_cone,
                            load_solution, forward, eqcode, eqvars, lines, NEQ, NVARS)
sys.setrecursionlimit(1000000)
CORE = {2071,4573,7123,7469,11854,13660,15299,16622,17726,21382,22093,25480,25539,28653,29437,31061,32894,32916,34517,34892}
gvids = {t: gates[definer[t]][2] for t in order}
base = load_solution('best/new_instance_partial_39013.json'); forward(base)
def count(v):
    ns={'__builtins__':{},'v':v}; return set(i for i in range(NEQ) if eval(eqcode[i],ns)!=0)
allvarcount=Counter()
for i in range(NEQ):
    for x in eqvars[i]: allvarcount[x]+=1
const_pinned=set()
_bp=re.compile(r'\(x_(\d+)\)\s*-\s*\(?-?\d{20,}'); _bp2=re.compile(r'-?\d{20,}\)?\s*-\s*\(x_(\d+)\)')
for i in range(NEQ):
    L=lines[i]
    if len(L)<4000:
        for m in _bp.finditer(L): const_pinned.add(int(m.group(1)))
        for m in _bp2.finditer(L): const_pinned.add(int(m.group(1)))
_,frS=backward_cone(35389); _,frT=backward_cone(6671)
STcone=set(frS)|set(frT); QUOT={30317,2936,5146}
def free_cone(r):
    seen=set(); st=[r]
    while st:
        u=st.pop()
        if u in seen: continue
        seen.add(u)
        for w in gvids.get(u,()):
            if w not in seen: st.append(w)
    return [u for u in seen if u in freeinp]
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

protect=(STcone|QUOT|const_pinned)-{31339,14853}
protect|={22152,33462}

def run(seed):
    rng=random.Random(seed)
    v=base[:]
    v[31339]=base[31339]+(base[16742]-base[24908]); forward(v)
    v[14853]=v[12186]; forward(v); set_quot(v)
    _kc={}
    def knob(gate):
        if gate in _kc: return _kc[gate]
        g0=v[gate]; cands=[]
        for w in free_cone(gate):
            if w in protect: continue
            old=v[w]; v[w]=old+1; forward(v); dd=v[gate]-g0; v[w]=old
            if dd in (1,-1): cands.append((allvarcount[w],w,dd))
        forward(v)
        if not cands: _kc[gate]=(None,None); return None,None
        cands.sort()
        # random among the lowest-multiplicity third (bias to clean but explore)
        k=max(1,len(cands)//3); pick=rng.choice(cands[:k]) if seed>0 else cands[0]
        _kc[gate]=(pick[1],pick[2]); return pick[1],pick[2]
    dirty=[i for i in range(NVARS) if i in freeinp and v[i]!=base[i]]
    if seed>0: rng.shuffle(dirty)
    queue=deque(dirty); assigned=set(dirty); aligned=set()
    cur=count(v); best=(len(cur),None); rounds=0
    while queue and rounds<40000:
        rounds+=1
        f=queue.popleft(); fval=v[f]
        ps=list(set().union(*[partners_in(i,f) for i in eqs_with(f)]) if eqs_with(f) else set())
        if seed>0: rng.shuffle(ps)
        for G in ps:
            if G in freeinp or G in aligned: continue
            if v[G]==fval: continue
            w,dd=knob(G)
            if w is None: continue
            old=v[w]; v[w]=old+dd*(fval-v[G]); forward(v); set_quot(v)
            if v[G]!=fval: v[w]=old; forward(v); set_quot(v); continue
            aligned.add(G); cur=count(v)
            if w not in assigned: assigned.add(w); queue.append(w)
            if len(cur)<best[0]: best=(len(cur),[x for x in v])
        if not cur: best=(0,[x for x in v]); break
    return best

globbest=(NEQ,None)
for seed in range(12):
    b=run(seed)
    tag=""
    if b[1] is not None and b[0]<globbest[0]:
        globbest=b; tag=" <== NEW BEST"
    print(f"seed {seed}: best {NEQ-b[0]}/{NEQ} ({b[0]} fail){tag}",flush=True)
    if globbest[0]==0: break
bv=globbest[1]
if bv is not None:
    cur=count(bv)
    print(f"\nGLOBAL BEST: {NEQ-len(cur)}/{NEQ} ({len(cur)} fail) core={len(cur&CORE)}: {sorted(cur)}")
    if NEQ-len(cur)>39022:
        json.dump({f"x_{i}":bv[i] for i in range(NVARS) if bv[i]!=0}, open(f'best_agentA_{NEQ-len(cur)}.json','w'))
        print(f"SAVED best_agentA_{NEQ-len(cur)}.json")
