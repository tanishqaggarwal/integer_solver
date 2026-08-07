#!/usr/bin/env python3
"""Relocate the chain break and score; plus free-dial optimisation."""
import sys,os,pickle,json,collections,itertools,time
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE)
from fwd import Engine,NV
E=Engine(); NR=len(E.res)
p=115792089237316195423570985008687907853269984665640564039457584007908834671663

def jac_at(v0,knobs):
    v=list(v0); F0=E.run(v); cols={}
    for f in knobs:
        old=v[f]; v[f]=old+1; F1=E.run(v); v[f]=old
        c={i:F1[i]-F0[i] for i in range(NR) if F1[i]!=F0[i]}
        if c: cols[f]=c
    return F0,cols

def components(cols,seed):
    rows=collections.defaultdict(dict)
    for f,c in cols.items():
        for i,val in c.items(): rows[i][f]=val
    colrows=collections.defaultdict(set)
    for i,r in rows.items():
        for f in r: colrows[f].add(i)
    seen=set(); comps=[]
    for s in seed:
        if s in seen: continue
        stack=[s]; comp=set()
        while stack:
            i=stack.pop()
            if i in seen: continue
            seen.add(i); comp.add(i)
            for f in rows.get(i,{}):
                for j in colrows[f]:
                    if j not in seen: stack.append(j)
        comps.append(sorted(comp))
    return rows,colrows,comps

def solve_comp(rows,comp,rhs,skip):
    """Solve rows in comp except `skip` exactly over Z; return delta or None."""
    active={i for i in comp if i!=skip}
    colr=collections.defaultdict(set)
    for i in active:
        for f in rows[i]: colr[f].add(i)
    delta={}
    changed=True
    while changed and active:
        changed=False
        for i in list(active):
            r=rows[i]
            unk=[f for f in r if f not in delta]
            if not unk:
                if sum(r[f]*delta[f] for f in r)==rhs[i]: active.discard(i); changed=True
                continue
            if len(unk)==1:
                f=unk[0]
                if any(j in active and j!=i for j in colr[f]): continue
                b=rhs[i]-sum(r[g]*delta[g] for g in r if g in delta)
                if b % r[f]: continue
                delta[f]=b//r[f]; active.discard(i); changed=True
    if active: return None
    return delta
