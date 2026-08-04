#!/usr/bin/env python3
"""Meta-iteration: run core-fix cascade; detect 'bad knobs' (clean frees the cascade moved
that break an equality-check whose partner gate has no knob); protect them; re-run. The
baseline satisfied ALL these equality-checks, so protecting bad knobs localizes the fix."""
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
QUOT={30317,2936,5146}
protect0=(set(frS)|set(frT)|QUOT|const_pinned|{22152,33462,16742,12186,24908,3558,29322,31339,14853})
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

def run(extra_protect):
    protect=protect0|extra_protect
    _kc={}
    def additive_knob(v,gate):
        if gate in _kc: return _kc[gate]
        g0=v[gate]; cands=[]
        for w in free_cone(gate):
            if w in protect: continue
            old=v[w]; v[w]=old+1; forward(v); dd=v[gate]-g0; v[w]=old
            if dd in (1,-1): cands.append((allvarcount[w],w,dd))
        forward(v); cands.sort()
        res=(cands[0][1],cands[0][2]) if cands else (None,None); _kc[gate]=res; return res
    v=base[:]
    v[31339]=base[31339]+(base[16742]-base[24908]); forward(v)
    v[14853]=v[12186]; forward(v); set_quot(v)
    dirty=[i for i in range(NVARS) if i in freeinp and v[i]!=base[i] and i not in protect]
    queue=deque(dirty); assigned=set(dirty); aligned=set(); cur=count(v); best=(len(cur),[x for x in v]); rounds=0
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
        if not cur: break
    v=best[1]; cur=count(v)
    # detect bad knobs: dirty clean frees w that appear in a failing eq as (w - gate)/(gate - w)
    # with baseline w==gate (equality-check the heal broke) and gate has no knob
    bad=set()
    for c in cur:
        for m in re.finditer(r'\(x_(\d+)\)\s*-\s*\(x_(\d+)\)',lines[c]):
            a,b=int(m.group(1)),int(m.group(2))
            for w,g in [(a,b),(b,a)]:
                if w in freeinp and w not in protect and v[w]!=base[w] and base[w]==base[g]:
                    bad.add(w)
    return v,cur,bad

extra=set(); best_overall=(len([]),None)
prev=None
for it in range(15):
    v,cur,bad=run(extra)
    print(f"meta {it}: {NEQ-len(cur)}/{NEQ} ({len(cur)} fail) core={len(cur&CORE)} extra_protect={len(extra)} new_bad={sorted(bad-extra)[:10]}",flush=True)
    if not cur:
        best_overall=(0,[x for x in v]); print("ALL SOLVED!"); break
    if best_overall[1] is None or len(cur)<len([i for i in range(NEQ) if False]) or (best_overall[0]==0):
        pass
    # track best
    if best_overall[1] is None or len(cur) < best_overall[0] if isinstance(best_overall[0],int) and best_overall[0]>0 else True:
        pass
    curbest = len(cur)
    if best_overall[1] is None or curbest < best_overall[0]:
        best_overall=(curbest,[x for x in v])
    newbad=bad-extra
    if not newbad:
        print("no new bad knobs; converged"); break
    extra|=newbad
bv=best_overall[1];
if bv is not None:
    cur=count(bv)
    print(f"BEST: {NEQ-len(cur)}/{NEQ} ({len(cur)} fail): {sorted(cur)}")
    if NEQ-len(cur)>39022:
        json.dump({f"x_{i}":bv[i] for i in range(NVARS) if bv[i]!=0}, open(f'best_agentA_{NEQ-len(cur)}.json','w'))
        print(f"SAVED best_agentA_{NEQ-len(cur)}.json")
