"""WR step 16: REALISE the uniform-plus-copy-knob deformation and measure it.

d = (1-p)*1 + sum_j t_j e_j  with t supported on wire coordinates chosen to zero
as many of the twelve a37694 rows as possible.  Integer solution found by HNF.
"""
import os, sys, collections, itertools, json, math
HERE = '/home/user/integer_solver/solve_lab/s10'
LAB = '/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
import wr_frame as W
import wr_rows as R
P = ad.P
WIRE, widx, rows, RE = R.WIRE, R.widx, R.rows, R.RE
N = len(WIRE)
rowsum = {e: sum(rows[e].values()) for e in RE}
BAD12 = [e for e in RE if rowsum[e]]

DET = dict(W.F3)
for u in WIRE:
    if u in L.definer:
        DET[u] = L.definer[u]
FW = W.Frame(DET)
base = L.load(os.path.join(LAB, 'best', 'new_instance_partial_39026.json'))
b2 = list(base); FW.fwd(b2)


def int_solve(A, b):
    """Integer solution t of A t = b (A: m x n ints).  Column HNF.  None if none."""
    m, n = len(A), len(A[0])
    M = [row[:] for row in A]
    U = [[1 if i == j else 0 for j in range(n)] for i in range(n)]
    piv = []
    r = 0
    for i in range(m):
        while True:
            nz = [c for c in range(n) if c not in piv and M[i][c]]
            if len(nz) <= 1:
                break
            nz.sort(key=lambda c: abs(M[i][c]))
            p0 = nz[0]
            for c in nz[1:]:
                q = M[i][c] // M[i][p0]
                if q:
                    for k in range(m):
                        M[k][c] -= q * M[k][p0]
                    for k in range(n):
                        U[k][c] -= q * U[k][p0]
        nz = [c for c in range(n) if c not in piv and M[i][c]]
        if nz:
            piv.append(nz[0])
    # now M is column-echelon on the pivot columns
    t = [0] * n
    bb = b[:]
    for i in range(m):
        pc = None
        for c in piv:
            if M[i][c]:
                # ensure this pivot is only used for its own row
                pc = c
                break
        if pc is None:
            if bb[i]:
                return None
            continue
        if bb[i] % M[i][pc]:
            return None
        q = bb[i] // M[i][pc]
        for k in range(m):
            bb[k] -= q * M[k][pc]
        for k in range(n):
            t[k] += q * U[k][pc]
    if any(bb):
        return None
    return t


def broken_rows(dvec):
    bad = []
    for e in RE:
        s = 0
        for j, c in rows[e].items():
            if dvec[j]:
                s += c * dvec[j]
        if s:
            bad.append(e)
    return bad


def measure(dvec, tag):
    v = list(b2)
    for j, u in enumerate(WIRE):
        v[u] = P + dvec[j]
    FW.fwd(v, rounds=10)
    av = L.all_atom_values(v)
    nz = [a for a in range(L.NA) if av[a]]
    fail = L.failing_eqs(av)
    print(f'{tag}: score={L.NEQ-len(fail)} failing={len(fail)} nonzero={len(nz)} '
          f'{sorted(nz)[:25]}', flush=True)
    print(f'   failing: {sorted(fail)}')
    return v, nz, fail


if __name__ == '__main__':
    cand = sorted(set().union(*[set(rows[e]) for e in BAD12]))
    print(f'candidate coordinates: {len(cand)} -> {[WIRE[j] for j in cand]}')
    C = 1 - P                              # wire value 1
    best = None
    for k in range(1, 7):
        loc = None
        for tgt in itertools.combinations(BAD12, k):
            A = [[rows[e].get(j, 0) for j in cand] for e in tgt]
            b = [-C * rowsum[e] for e in tgt]
            t = int_solve(A, b)
            if t is None:
                continue
            d = [C] * N
            for i, j in enumerate(cand):
                d[j] += t[i]
            bad = broken_rows(d)
            if loc is None or len(bad) < len(loc[0]):
                loc = (bad, d, tgt)
        if loc is None:
            print(f'k={k}: no integer solution')
            continue
        print(f'k={k}: {len(loc[0])} identity rows broken, targets {loc[2]}')
        print(f'      broken {loc[0]}')
        if best is None or len(loc[0]) < len(best[0]):
            best = loc
    if best:
        bad, d, tgt = best
        mx = max(abs(x) for x in d)
        print(f'\nBEST: {len(bad)} identity rows; max |d| digits {len(str(mx))}')
        v, nz, fail = measure(d, 'uniform+knobs')
        T.save(v, os.path.join(HERE, 'wr_knobs.json'))
