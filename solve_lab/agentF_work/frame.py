#!/usr/bin/env python3
"""Frame analysis: given a set S of atoms allowed to be nonzero, compute the
integrally reachable residual lattice and the max number of satisfiable equations."""
import sys,os,json,itertools,collections,time
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE)
from intsolve import solve_int

def analyze(F,S,maxrows=22,verbose=False):
    Sset=set(S); S=sorted(S)
    E=sorted({e for j in S for e,k in F.atom2eq[j]})
    if len(E)>maxrows: return None
    priv=[i for i in range(len(F.v)) if F.var2atoms[i] and set(F.var2atoms[i])<=Sset]
    if not priv: return None
    base=[F.av[j] for j in S]
    G=[]
    for i in priv:
        old=F.v[i]; F.setvar(i,old+1)
        col=[F.av[j]-base[k] for k,j in enumerate(S)]
        F.setvar(i,old)
        G.append(col)
    # G is list of columns (len(priv) columns, each len(S))
    M=[[dict((jj,kk) for kk,jj in F.eqrows[e]).get(j,0) for j in S] for e in E]
    nq=len(priv)
    A_full=[[sum(M[r][t]*G[c][t] for t in range(len(S))) for c in range(nq)] for r in range(len(E))]
    b_full=[-sum(M[r][t]*base[t] for t in range(len(S))) for r in range(len(E))]
    best=(0,None,None)
    n=len(E)
    for size in range(n,0,-1):
        found=False
        for T in itertools.combinations(range(n),size):
            x=solve_int([A_full[r] for r in T],[b_full[r] for r in T])
            if x is not None:
                best=(size,[E[r] for r in T],dict(zip(priv,x))); found=True; break
        if found: break
    return dict(S=S,E=E,priv=priv,nsat=best[0],rows=best[1],sol=best[2],nfail=len(E)-best[0])
