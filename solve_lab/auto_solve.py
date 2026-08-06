#!/usr/bin/env python3
"""Auto-expanding coupled integer solver for the twist.

Iteratively: (1) grow the unknown set with vars that appear only-linearly-given-UNK
(never multiplied by another UNK var), (2) put EVERY atom touching a UNK var into the
system, (3) solve A x = b over Z (SNF), (4) scan all atoms; add any still-nonzero atom's
frontier and repeat. Product-partners stay at base so the system is linear. Verify with
the checker."""
import json, sys
from collections import defaultdict
from propagate import load_atoms, atom_vars, NVARS

def snf_solve(A, b):
    m=len(A); n=len(A[0])
    M=[r[:] for r in A]; bb=b[:]
    V=[[1 if i==j else 0 for j in range(n)] for i in range(n)]
    for piv in range(min(m,n)):
        while True:
            best=None; found=None
            for i in range(piv,m):
                row=M[i]
                for j in range(piv,n):
                    if row[j]!=0 and (best is None or abs(row[j])<best):
                        best=abs(row[j]); found=(i,j)
            if found is None: break
            i0,j0=found
            if i0!=piv: M[piv],M[i0]=M[i0],M[piv]; bb[piv],bb[i0]=bb[i0],bb[piv]
            if j0!=piv:
                for row in M: row[piv],row[j0]=row[j0],row[piv]
                for row in V: row[piv],row[j0]=row[j0],row[piv]
            p=M[piv][piv]; done=True
            for i in range(m):
                if i!=piv and M[i][piv]!=0:
                    q=M[i][piv]//p
                    Mi=M[i]; Mp=M[piv]
                    for j in range(n): Mi[j]-=q*Mp[j]
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
    for i in range(min(m,n)):
        d=M[i][i]
        if d==0: continue
        if bb[i]%d!=0: return None
        y[i]=bb[i]//d
    for i in range(m):
        if all(v==0 for v in M[i]) and bb[i]!=0: return None
    return [sum(V[i][j]*y[j] for j in range(n)) for i in range(n)]

A=load_atoms()
base={int(k[2:]):x for k,x in json.load(open('rebuilt_partial.json')).items()}
va=defaultdict(list)
for i,poly in enumerate(A):
    for x in atom_vars(poly): va[x].append(i)
def ev(poly,val):
    s=0
    for m_,c in poly.items():
        t=c
        for z in m_: t*=val.get(z,0)
        s+=t
    return s

# seed unknowns: the free product partners + spine flow vars
UNK=set([14257,18956,7497,32237,37892,10156,25538,22820,
         24468,13682,34243,13913,38045,14393,11436,32989,9274,29237,23134])
cur=dict(base)

def safe_to_add(z, UNK):
    # z must never be multiplied by another UNK var (keep system linear)
    if abs(base.get(z,0))>=10**6: return False   # already-huge non-flow var: don't touch
    for ai in va[z]:
        for m_ in A[ai]:
            if z in m_ and len(m_)>=2:
                for w in m_:
                    if w!=z and w in UNK: return False
    return True

for it in range(25):
    # system = all atoms touching any UNK var
    SYS=set()
    for z in UNK: SYS.update(va[z])
    SYS=sorted(SYS)
    # grow UNK with linear-only vars from these atoms (fixpoint-ish, one pass)
    add=set()
    for ai in SYS:
        for z in atom_vars(A[ai]):
            if z not in UNK and safe_to_add(z,UNK): add.add(z)
    if add:
        UNK|=add
        continue
    idx={v:i for i,v in enumerate(sorted(UNK))}; U=sorted(UNK); n=len(U)
    rows=[]; bs=[]
    for ai in SYS:
        coef=[0]*n; const=0
        for m_,c in A[ai].items():
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
        print(f"it{it}: no integer solution ({len(SYS)} atoms, {n} unk)"); break
    cur=dict(base)
    for v in U: cur[v]=x[idx[v]]
    nz=[i for i in range(len(A)) if ev(A[i],cur)!=0]
    print(f"it{it}: {len(SYS)} atoms, {n} unk -> global nonzero atoms {len(nz)}", flush=True)
    if not nz:
        json.dump({f"x_{i}":cur.get(i,0) for i in range(NVARS)}, open('auto_solved.json','w'))
        print("*** ALL ATOMS ZERO -> wrote auto_solved.json ***"); break
    # expand: add frontier of nonzero atoms
    newu=set()
    for ai in nz:
        for z in atom_vars(A[ai]):
            if z not in UNK and safe_to_add(z,UNK): newu.add(z)
    if not newu:
        print(f"  cannot expand further; {len(nz)} atoms stuck: {nz[:12]}")
        json.dump({f"x_{i}":cur.get(i,0) for i in range(NVARS)}, open('auto_solved.json','w'))
        break
    UNK|=newu
json.dump({f"x_{i}":cur.get(i,0) for i in range(NVARS)}, open('auto_solved.json','w'))
