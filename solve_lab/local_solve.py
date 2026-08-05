#!/usr/bin/env python3
"""Solve the residual twist as a coupled INTEGER-linear system (Smith normal form).

Two product-slack chains (multipliers = V0-pinned wire) couple through shared verifier
checks. Build the local system, hold vars outside the unknown set at rebuilt_partial
(wire = V0, small product-partners fixed => every monomial linear), and solve A x = b
over Z via SNF. Verify system atoms exactly, then run the checker."""
import json, sys
from collections import defaultdict
from propagate import load_atoms, atom_vars, NVARS

def snf_solve(A, b):
    """Integer solution of A x = b (A: m lists of n ints). Returns x (n ints) or None."""
    m = len(A); n = len(A[0])
    M = [row[:] for row in A]; bb = b[:]
    V = [[1 if i==j else 0 for j in range(n)] for i in range(n)]
    r = min(m, n)
    for piv in range(r):
        while True:
            best=None; found=None
            for i in range(piv, m):
                for j in range(piv, n):
                    if M[i][j]!=0 and (best is None or abs(M[i][j])<best):
                        best=abs(M[i][j]); found=(i,j)
            if found is None: break
            i0,j0=found
            if i0!=piv:
                M[piv],M[i0]=M[i0],M[piv]; bb[piv],bb[i0]=bb[i0],bb[piv]
            if j0!=piv:
                for row in M: row[piv],row[j0]=row[j0],row[piv]
                for row in V: row[piv],row[j0]=row[j0],row[piv]
            p=M[piv][piv]; done=True
            for i in range(m):
                if i!=piv and M[i][piv]!=0:
                    q=M[i][piv]//p
                    for j in range(n): M[i][j]-=q*M[piv][j]
                    bb[i]-=q*bb[piv]
                    if M[i][piv]!=0: done=False
            for j in range(n):
                if j!=piv and M[piv][j]!=0:
                    q=M[piv][j]//p
                    for i in range(m): M[i][j]-=q*M[i][piv]
                    for i in range(n): V[i][j]-=q*V[i][piv]
                    if M[piv][j]!=0: done=False
            if done: break
    y=[0]*n
    rank=0
    for i in range(min(m,n)):
        d=M[i][i]
        if d==0: continue
        rank+=1
        if bb[i]%d!=0: return None
        y[i]=bb[i]//d
    for i in range(m):
        if all(M[i][j]==0 for j in range(n)) and bb[i]!=0: return None
    x=[sum(V[i][j]*y[j] for j in range(n)) for i in range(n)]
    return x

A = load_atoms()
base = {int(k[2:]): x for k, x in json.load(open('rebuilt_partial.json')).items()}
SYSTEM = [602,601,29375,29374,17686, 1465,1464,29373,29372,17681,
          40907,44255,46225,45828, 21603,21602,16554]
UNK = {14257,18956,7497,32237,37892,10156,25538,22820,
       24468,13682,34243,13913,38045,14393,11436,32989,
       9274,29237,23134}
# add check vars that appear ONLY in linear monomials (safe to adjust; product-partners stay 0)
for ai in (40907,44255,46225,45828,21602):
    lin=set(); prod=set()
    for m_ in A[ai]:
        if len(m_)==1: lin.add(m_[0])
        elif len(m_)>=2:
            for z in m_: prod.add(z)
    for z in lin-prod:
        if abs(base.get(z,0))<10**6: UNK.add(z)
UNK = sorted(UNK); idx={v:i for i,v in enumerate(UNK)}; n=len(UNK)
print(f"{len(SYSTEM)} atoms, {n} unknowns", flush=True)

rows=[]; bs=[]
for ai in SYSTEM:
    coef=[0]*n; const=0
    for m_, c in A[ai].items():
        us=[x for x in m_ if x in idx]
        if not us:
            t=c
            for x in m_: t*=base.get(x,0)
            const+=t
        else:
            keep=us[0]; t=c
            for x in m_:
                if x==keep: continue
                t*=base.get(x,0)
            coef[idx[keep]]+=t
    rows.append(coef); bs.append(-const)

x=snf_solve(rows,bs)
if x is None:
    print("NO integer solution to the linearized system"); sys.exit(0)
cand=dict(base)
for v in UNK: cand[v]=x[idx[v]]
def ev(poly,val):
    s=0
    for m_,c in poly.items():
        t=c
        for z in m_: t*=val.get(z,0)
        s+=t
    return s
bad=[ai for ai in SYSTEM if ev(A[ai],cand)!=0]
print(f"system atoms still nonzero: {bad}", flush=True)
json.dump({f"x_{i}":cand.get(i,0) for i in range(NVARS)}, open('local_solved.json','w'))
print("wrote local_solved.json", flush=True)
