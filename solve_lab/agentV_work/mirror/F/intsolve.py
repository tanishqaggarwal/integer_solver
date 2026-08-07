#!/usr/bin/env python3
"""Exact integer solve A x = b (small dense) via column HNF."""
from fractions import Fraction

def col_hnf(A):
    """A: list of rows (lists). Returns (H,U) with A@U=H, U unimodular, H in column echelon."""
    m=len(A); n=len(A[0])
    H=[row[:] for row in A]
    U=[[1 if i==j else 0 for j in range(n)] for i in range(n)]
    def colop(j,k,q):  # col j -= q*col k
        for i in range(m): H[i][j]-=q*H[i][k]
        for i in range(n): U[i][j]-=q*U[i][k]
    def swap(j,k):
        for i in range(m): H[i][j],H[i][k]=H[i][k],H[i][j]
        for i in range(n): U[i][j],U[i][k]=U[i][k],U[i][j]
    piv=0; pivots=[]
    for r in range(m):
        # find nonzero in row r among cols piv..n-1, reduce via gcd
        while True:
            nz=[j for j in range(piv,n) if H[r][j]!=0]
            if not nz: break
            if len(nz)==1: break
            nz.sort(key=lambda j: abs(H[r][j]))
            j0=nz[0]
            done=True
            for j in nz[1:]:
                q=H[r][j]//H[r][j0]
                if q: colop(j,j0,q); done=False
            if all(H[r][j]==0 for j in range(piv,n) if j!=j0): break
            if done: break
        nz=[j for j in range(piv,n) if H[r][j]!=0]
        if not nz: pivots.append(None); continue
        j0=nz[0]
        if j0!=piv: swap(piv,j0)
        pivots.append(piv); piv+=1
    return H,U,pivots

def solve_int(A,b):
    """Return integer x with A x = b, or None."""
    m=len(A); n=len(A[0])
    H,U,pivots=col_hnf(A)
    y=[0]*n
    bb=b[:]
    for r in range(m):
        j=pivots[r]
        if j is None:
            if bb[r]!=0: return None
            continue
        if H[r][j]==0:
            if bb[r]!=0: return None
            continue
        if bb[r] % H[r][j]: return None
        y[j]=bb[r]//H[r][j]
        for i in range(m): bb[i]-=y[j]*H[i][j]
    if any(v!=0 for v in bb): return None
    x=[sum(U[i][j]*y[j] for j in range(n)) for i in range(n)]
    # verify
    for i in range(m):
        if sum(A[i][j]*x[j] for j in range(n))!=b[i]: return None
    return x
