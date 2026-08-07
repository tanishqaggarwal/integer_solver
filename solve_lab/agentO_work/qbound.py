"""Rational upper bound on maxsat for every single-atom growth of the witness region.

Integer solvability implies rational solvability, so max over Q is an upper bound on max
over Z, hence |E| - maxsat_Q is a LOWER bound on the number of failing equations.  Rank
tests are cheap, so this closes the single-atom growth question rigorously.
"""
import sys, json, time, itertools
from fractions import Fraction
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentO_work')
import regiongrow as G, harness as H

OD = '/home/user/integer_solver/solve_lab/agentO_work'
LOG = open(OD + '/runs/qbound.log', 'w', buffering=1)


def say(*a):
    print(*a, file=LOG)


def rank(rows):
    R = [list(map(Fraction, r)) for r in rows]
    n = len(R)
    m = len(R[0]) if n else 0
    rk = 0
    for c in range(m):
        p = None
        for i in range(rk, n):
            if R[i][c]:
                p = i
                break
        if p is None:
            continue
        R[rk], R[p] = R[p], R[rk]
        pv = R[rk][c]
        for i in range(n):
            if i != rk and R[i][c]:
                f = R[i][c] / pv
                for k in range(c, m):
                    R[i][k] -= f * R[rk][k]
        rk += 1
    return rk


def qsat(S, A, B, order):
    rows = [[A[e].get(u, 0) for u in order] for e in S]
    aug = [rows[i] + [B[S[i]]] for i in range(len(S))]
    return rank(rows) == rank(aug)


def qmaxsat(Eqs, A, B, order):
    Eqs = list(Eqs)
    idx = {e: i for i, e in enumerate(Eqs)}
    cur = {(e,) for e in Eqs if qsat([e], A, B, order)}
    best = 1 if cur else 0
    bestS = list(next(iter(cur))) if cur else []
    while cur:
        nxt = set()
        for S in cur:
            last = idx[S[-1]]
            for e in Eqs[last + 1:]:
                T = S + (e,)
                if not all(T[:j] + T[j + 1:] in cur for j in range(len(T))):
                    continue
                if qsat(list(T), A, B, order):
                    nxt.add(T)
        if not nxt:
            break
        cur = nxt
        T = next(iter(cur))
        best = len(T)
        bestS = list(T)
    return best, bestS


cands = json.load(open(OD + '/growcand.json'))
say('candidate atoms that free a private var: %d' % len(cands))
rows_out = []
for c in cands:
    a = c['atom']
    R = G.R0 + [a]
    P = G.private_vars(R)
    m = G.build_model(R, P, G.V0)
    if m is None:
        say('  a%d: NONLINEAR' % a)
        continue
    const, cols = m
    Eqs, rows = G.eq_system(R, P, const, cols)
    A = {e: rows[e][0] for e in Eqs}
    B = {e: rows[e][1] for e in Eqs}
    t0 = time.time()
    k, S = qmaxsat(Eqs, A, B, P)
    lb = len(Eqs) - k
    rows_out.append((lb, a, len(Eqs), k, len(P)))
    say('  a%d: P=%d |E|=%d  maxsat_Q=%d  cost >= %d  (%.0fs)%s'
        % (a, len(P), len(Eqs), k, lb, time.time() - t0, '   <<< could beat 7' if lb < 7 else ''))
rows_out.sort()
say('sorted by lower bound: %s' % (rows_out[:15],))
json.dump(rows_out, open(OD + '/qbound.json', 'w'))
say('DONE')
