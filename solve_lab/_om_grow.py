#!/usr/bin/env python3
"""Grow the atom/equation closure around the two bad atoms; track rank & whether
any null vector can have t_alpha != 0  (mod a big prime, sparse elimination)."""
import pickle, sys, random
from collections import defaultdict

Pm = (1 << 61) - 1  # prime modulus for rank computation

D = pickle.load(open('_om_parsed2.pkl', 'rb'))
eqatoms = D['eqatoms']
ainc = defaultdict(list)
for e, d in enumerate(eqatoms):
    for k in d: ainc[k].append(e)
F = [2554, 6816, 8124, 8680, 9421, 12231, 12270, 12350, 14584, 22044, 29125]
ALPHA = '((x7068-x2099)-(7376877*x642))'
BETA = '((x4432-x19964)-x28730)'

def sparse_rank(rows, cols_order):
    """rows: list of dict col->coeff (ints).  Returns (rank, pivot_cols).
    Sparse gaussian elimination mod Pm with simple Markowitz-ish ordering."""
    rows = [{c: v % Pm for c, v in r.items() if v % Pm} for r in rows]
    rows = [r for r in rows if r]
    pivots = {}   # col -> row dict (normalized, pivot coeff 1)
    order = {c: i for i, c in enumerate(cols_order)}
    work = sorted(rows, key=len)
    import heapq
    heap = [(len(r), i) for i, r in enumerate(work)]
    heapq.heapify(heap)
    alive = list(work)
    n = 0
    while heap:
        ln, i = heapq.heappop(heap)
        r = alive[i]
        if r is None or len(r) != ln:
            if r is not None: heapq.heappush(heap, (len(r), i))
            continue
        # reduce by existing pivots
        changed = True
        while changed:
            changed = False
            for c in list(r.keys()):
                if c in pivots:
                    f = r[c]
                    pr = pivots[c]
                    for cc, vv in pr.items():
                        nv = (r.get(cc, 0) - f * vv) % Pm
                        if nv: r[cc] = nv
                        elif cc in r: del r[cc]
                    changed = True
                    break
        if not r:
            alive[i] = None
            continue
        # pick pivot = least-frequent column
        pc = min(r.keys(), key=lambda c: (order.get(c, 1 << 30), c))
        inv = pow(r[pc], Pm - 2, Pm)
        pr = {c: (v * inv) % Pm for c, v in r.items()}
        pivots[pc] = pr
        alive[i] = None
        n += 1
    return n, set(pivots.keys())

def closure(seed_atoms, shells):
    A = set(seed_atoms)
    for s in range(shells):
        E = set(F)
        for k in A: E |= set(ainc[k])
        newA = set(A)
        for e in E: newA |= set(eqatoms[e])
        if newA == A: break
        A = newA
        E2 = set(F)
        for k in A: E2 |= set(ainc[k])
        print('shell %d: |A|=%d |E|=%d' % (s + 1, len(A), len(E2)))
    E = set(F)
    for k in A: E |= set(ainc[k])
    return sorted(A), sorted(E)

if __name__ == '__main__':
    ns = int(sys.argv[1]) if len(sys.argv) > 1 else 2
    A0 = set()
    for e in F: A0 |= set(eqatoms[e])
    A, E = closure(A0, ns)
    print('final |A|=%d |E|=%d' % (len(A), len(E)))
    colfreq = defaultdict(int)
    for e in E:
        for k in eqatoms[e]:
            if k in set(A): colfreq[k] += 1
    Aset = set(A)
    rows = [{k: c for k, c in eqatoms[e].items() if k in Aset} for e in E]
    order = sorted(A, key=lambda k: colfreq[k])
    rk, piv = sparse_rank(rows, order)
    print('rank=%d  nullity=%d' % (rk, len(A) - rk))
    # Can t_alpha be nonzero?  add row  t_alpha = 0 and see if rank grows
    rk2, _ = sparse_rank(rows + [{ALPHA: 1}], order)
    rk3, _ = sparse_rank(rows + [{BETA: 1}], order)
    print('rank with t_alpha=0 forced: %d  (grew=%s -> alpha CAN be nonzero: %s)' % (rk2, rk2 > rk, rk2 > rk))
    print('rank with t_beta =0 forced: %d  (grew=%s -> beta  CAN be nonzero: %s)' % (rk3, rk3 > rk, rk3 > rk))
