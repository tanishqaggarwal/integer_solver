#!/usr/bin/env python3
"""Exact integer solve of the sparse residual system J*delta = -F0 by peeling."""
import sys,os,pickle,collections,json,time
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE)

def solve_peel(rows, rhs):
    """rows: {i:{f:c}}, rhs: {i:b}. Returns (delta, unsolved_rows, notes)."""
    rows={i:dict(r) for i,r in rows.items()}
    colrows=collections.defaultdict(set)
    for i,r in rows.items():
        for f in r: colrows[f].add(i)
    delta={}
    active=set(rows)
    changed=True
    while changed:
        changed=False
        for i in list(active):
            r=rows[i]
            unk=[f for f in r if f not in delta]
            if len(unk)==0:
                s=sum(r[f]*delta[f] for f in r)
                if s==rhs[i]: active.discard(i); changed=True
                continue
            if len(unk)==1:
                f=unk[0]
                # is f exclusive to active rows containing it?
                others=[j for j in colrows[f] if j in active and j!=i]
                if others: continue
                b=rhs[i]-sum(r[g]*delta[g] for g in r if g in delta)
                c=r[f]
                if b % c != 0: continue
                delta[f]=b//c
                active.discard(i); changed=True
    return delta, sorted(active)
