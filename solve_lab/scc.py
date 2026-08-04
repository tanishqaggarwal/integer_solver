#!/usr/bin/env python3
"""Build canonical gate DAG from gates.jsonl. Find SCCs (cycles) = the hard core.
Each target may have multiple defs (join points); we consider ALL of them for cycle detection
(a var depends on the union of inputs across its definitions)."""
import json, sys
from collections import defaultdict
sys.setrecursionlimit(1000000)

gates=[]
with open('atoms/gates.jsonl') as f:
    for line in f:
        d=json.loads(line); gates.append((d['t'], d['vids']))
# target -> set of input vars (union over all defs)
deps=defaultdict(set)
targets=set()
for t,vids in gates:
    targets.add(t)
    for v in vids: deps[t].add(v)
allvars=set(deps)|set(v for vs in deps.values() for v in vs)
print(f"targets={len(targets)} vars touched={len(allvars)}")

# Tarjan SCC (iterative)
index={}; low={}; onstk={}; stk=[]; idx=[0]; sccs=[]
def strongconnect(v):
    work=[(v,0)]
    while work:
        node,pi=work[-1]
        if pi==0:
            index[node]=low[node]=idx[0]; idx[0]+=1; stk.append(node); onstk[node]=True
        recurse=False
        neigh=list(deps.get(node,()))
        i=pi
        while i<len(neigh):
            w=neigh[i]
            if w not in index:
                work.append((w,0)); work[-2]=(node,i+1); recurse=True; break
            elif onstk.get(w): low[node]=min(low[node],index[w])
            i+=1
        if recurse: continue
        for w in neigh:
            if w in low and onstk.get(w): low[node]=min(low[node],low[w])
        if low[node]==index[node]:
            comp=[]
            while True:
                w=stk.pop(); onstk[w]=False; comp.append(w)
                if w==node: break
            sccs.append(comp)
        work.pop()
        if work:
            par=work[-1][0]; low[par]=min(low[par],low[node])
for v in list(allvars):
    if v not in index: strongconnect(v)
nontrivial=[c for c in sccs if len(c)>1]
# also self-loops (var depends on itself directly)
selfloop=[t for t in deps if t in deps[t]]
print(f"total SCCs={len(sccs)} nontrivial (cycles)={len(nontrivial)} self-loops={len(selfloop)}")
nontrivial.sort(key=len, reverse=True)
print("largest cycles:", [len(c) for c in nontrivial[:10]])
if nontrivial:
    big=nontrivial[0]
    print(f"biggest cycle ({len(big)} vars):", sorted(big)[:40])
json.dump({'nontrivial':[sorted(c) for c in nontrivial], 'selfloop':selfloop}, open('sccs.json','w'))
