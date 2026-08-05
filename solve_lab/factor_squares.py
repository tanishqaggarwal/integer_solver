#!/usr/bin/env python3
"""CREATIVE: does the verifier square root Q (deg-2) FACTOR into two linear forms
L1*L2?  A quadratic polynomial factors iff its (n+1)x(n+1) symmetric coefficient
matrix has rank <= 2.  If Q40782 = L1*L2, then a40782 = L1^2 * L2^2 and the verifier
square is satisfied by the LINEAR disjunction (L1=0 OR L2=0) -- collapsing the hard
quadratic into linear cases.  Compute the rank of the quadratic-form matrix over
GF(P) for a40782, a39550, and every perfect-square atom; report which factor."""
import json, time
from collections import defaultdict
from propagate import load_atoms, atom_vars
from check_square import try_sqrt
from modp import P, inv

def qform_rank(Q):
    """Q: dict monomial(tuple len<=2)->coef. Build symmetric matrix over GF(P) with
    index 0 = constant '1'. Return (rank, varlist)."""
    vs = sorted(atom_vars(Q))
    idx = {v: i+1 for i, v in enumerate(vs)}
    n = len(vs)+1
    M = [[0]*n for _ in range(n)]
    half = inv(2)
    for m, c in Q.items():
        c %= P
        if len(m) == 0:
            M[0][0] = (M[0][0]+c) % P
        elif len(m) == 1:
            i = idx[m[0]]
            M[0][i] = (M[0][i]+c*half) % P; M[i][0] = M[0][i]
        else:
            a, b = m
            if a == b:
                i = idx[a]; M[i][i] = (M[i][i]+c) % P
            else:
                i, j = idx[a], idx[b]
                M[i][j] = (M[i][j]+c*half) % P; M[j][i] = M[i][j]
    # gaussian rank over GF(P)
    r = 0; rows = [row[:] for row in M]
    col = 0; nr = n
    pr = 0
    for c in range(n):
        piv = None
        for rr in range(pr, nr):
            if rows[rr][c] % P: piv = rr; break
        if piv is None: continue
        rows[pr], rows[piv] = rows[piv], rows[pr]
        ipv = inv(rows[pr][c])
        for rr in range(nr):
            if rr != pr and rows[rr][c] % P:
                f = (rows[rr][c]*ipv) % P
                for cc in range(n):
                    rows[rr][cc] = (rows[rr][cc]-f*rows[pr][cc]) % P
        pr += 1
        if pr == nr: break
    return pr, vs

def main():
    t0 = time.time()
    A = load_atoms()
    # find all perfect-square atoms (deg-4 that are Q^2) + known ones
    print("scanning for perfect-square atoms...", flush=True)
    squares = []
    for a, poly in enumerate(A):
        if not poly: continue
        deg = max(len(m) for m in poly)
        if deg == 4:
            Q = try_sqrt(poly)
            if Q: squares.append((a, Q))
    print(f"found {len(squares)} deg-4 perfect squares ({time.time()-t0:.0f}s)", flush=True)
    lowrank = []
    for a, Q in squares:
        r, vs = qform_rank(Q)
        if r <= 2:
            lowrank.append((a, r, len(vs)))
    print(f"\nverifier squares whose root Q FACTORS into linear forms (rank<=2):", flush=True)
    for a, r, nv in lowrank:
        print(f"  a{a}: Q rank={r} ({nv} vars)  -> a{a}=0 is a LINEAR disjunction", flush=True)
    print(f"\n{len(lowrank)}/{len(squares)} verifier squares factor into linear forms", flush=True)
    # detail for the twist square a40782 and a39550
    for a in (40782, 39550):
        Q = try_sqrt(A[a]); r, vs = qform_rank(Q)
        print(f"  a{a}: Q quadratic-form rank = {r} (vars {len(vs)})  {'FACTORS' if r<=2 else 'irreducible quadric'}", flush=True)
    print(f"done ({time.time()-t0:.0f}s)", flush=True)

if __name__ == '__main__':
    main()
