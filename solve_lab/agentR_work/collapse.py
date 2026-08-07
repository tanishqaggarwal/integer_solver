#!/usr/bin/env python3
"""DOES THE GAP ACTUALLY COST ANYTHING?
tradeoff.py assumed the selectors BETWEEN the two relaxed ones are arbitrary booleans, so each
intervening stage applies the chord law and doubles the degree in t1.  But I choose those
selectors.  If every selector strictly between i and j (and after j) is 0, the mux
acc' = acc + b*(S-acc) is the IDENTITY there, acc does not move, and the degree does NOT grow.
Then acc_root = acc_j = A + t1*(S_i-A) + t2*(S_j-acc_i): 2 unknowns, 2 coordinates of T,
degree ~3 elimination -> trivially solvable, INDEPENDENT of the gap.
Tested on the scaled sibling of identical shape (small prime, brute force over t1)."""
import sys
import sibling

def test(m, i, j, verbose=True):
    d = sibling.instance(m)
    p, n, lad, T = d['p'], d['n'], d['lad'], d['T']
    def chord(A, B):
        ax, ay = A; bx, by = B
        if (ax - bx) % p == 0: return None
        l = (by - ay) * pow(bx - ax, p - 2, p) % p
        sx = (l * l - ax - bx) % p
        return (sx, (l * (ax - sx) - ay) % p)
    def mux(A, Sp, t):
        return ((A[0] + t * (Sp[0] - A[0])) % p, (A[1] + t * (Sp[1] - A[1])) % p)
    A = lad[0]                                  # accumulator seed, selectors 0..i-1 all OFF
    Si = chord(A, lad[i])
    if Si is None: return None
    hits = []
    for t1 in range(p):
        acc = mux(A, Si, t1)                    # stages between i and j: selectors 0 -> identity
        Sj = chord(acc, lad[j])
        if Sj is None: continue
        dx = (Sj[0] - acc[0]) % p; dy = (Sj[1] - acc[1]) % p
        if dx == 0: continue
        t2 = (T[0] - acc[0]) * pow(dx, p - 2, p) % p
        if (acc[1] + t2 * dy - T[1]) % p == 0:
            hits.append((t1, t2))
    if verbose:
        print('m=%-3d p=%-7d gap=%-4d relaxed leaves (%d,%d): %d solution(s) %s'
              % (m, p, j - i, i, j, len(hits), hits[:3]), flush=True)
    return hits

print('Question: with intervening selectors 0, does a (t1,t2) reaching T exist REGARDLESS of gap?')
for m in (8, 10, 12, 14, 16):
    d = sibling.instance(m); n = d['n']
    for (i, j) in ((1, 2), (1, n - 1), (2, n - 1)):
        if j <= i or j >= n: continue
        test(m, i, j)
