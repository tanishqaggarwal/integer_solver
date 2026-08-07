#!/usr/bin/env python3
"""Generic residual-atom repair: Gauss-Seidel over free inputs."""
import sys,os,pickle,collections,json,time,random
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE)
from fwd import Engine,NV
from circ2 import vars_of

E=Engine()
sup=pickle.load(open(os.path.join(HERE,'supp.pkl'),'rb'))
free=set(E.free)
# atom supports
asup=[]
for a in E.res:
    s=set()
    for u in vars_of(E.atoms[a]): s|=set(sup[str(u)])
    asup.append(frozenset(s))
reach=collections.defaultdict(set)
for i,s in enumerate(asup):
    for f in s: reach[f].add(i)
NR=len(E.res)

def ev(v):
    r=E.run(v); return r

def atom_val(v,i):
    # need full forward each time (cheap: 16ms)
    r=E.run(v); return r[i]

def full(v):
    r=E.run(v); bad=E.score(r); return r,bad

def nz_atoms(r): return [i for i in range(NR) if r[i]]

def try_fix(v,i,verbose=False):
    """Try to zero residual atom i by changing one free input in its support."""
    r=E.run(v); base=r[i]
    if base==0: return None
    cands=sorted(asup[i], key=lambda f: len(reach[f]))
    for f in cands:
        old=v[f]
        v[f]=old+1; r1=E.run(v); a1=r1[i]
        v[f]=old+2; r2=E.run(v); a2=r2[i]
        v[f]=old
        d1=a1-base; d2=a2-a1
        if d1==0 or d1!=d2: continue     # not affine in f
        if base % d1 != 0: continue
        step=-base//d1
        v[f]=old+step; r3=E.run(v)
        if r3[i]==0:
            return f,old,old+step
        v[f]=old
    return None
