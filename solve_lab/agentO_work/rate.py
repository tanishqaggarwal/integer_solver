"""STEP 1 — structure of the four blocking conditions at the witness configuration.

Compute the unique rational solution of the 13-equation region system exactly, read off the
denominators, and test whether the four blocking congruences are independent or coupled.
This decides whether a configuration scan is a search or a measurement.
"""
import sys, json, math
from fractions import Fraction
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentO_work')
import regiongrow as G, harness as H

OD = '/home/user/integer_solver/solve_lab/agentO_work'
LOG = open(OD + '/runs/rate.log', 'w', buffering=1)
P = 115792089237316195423570985008687907853269984665640564039457584007908834671663


def say(*a):
    print(*a, file=LOG)


R = G.R0 + [23618]
Pv = G.private_vars(R)
const, cols = G.build_model(R, Pv, G.V0)
Eqs, rows = G.eq_system(R, Pv, const, cols)
A = [[Fraction(rows[e][0].get(u, 0)) for u in Pv] for e in Eqs]
b = [Fraction(rows[e][1]) for e in Eqs]
say('unknowns:', Pv)
say('equations:', Eqs)

# exact rational solve
M = [A[i][:] + [b[i]] for i in range(len(A))]
n, m = len(M), len(Pv)
piv = []
r = 0
for c in range(m):
    p_ = next((i for i in range(r, n) if M[i][c]), None)
    if p_ is None:
        continue
    M[r], M[p_] = M[p_], M[r]
    pv = M[r][c]
    M[r] = [x / pv for x in M[r]]
    for i in range(n):
        if i != r and M[i][c]:
            f = M[i][c]
            M[i] = [M[i][k] - f * M[r][k] for k in range(m + 1)]
    piv.append(c)
    r += 1
say('rank over Q =', r, '(unknowns %d)' % m)
incons = [i for i in range(r, n) if M[i][m] != 0]
say('inconsistent dependent rows:', incons, '(empty => unique rational solution satisfies ALL)')

z = {}
for i, c in enumerate(piv):
    z[Pv[c]] = M[i][m]
say('\nunique rational solution — denominators:')
blocked = []
for u in Pv:
    q = z[u]
    d = q.denominator
    fac = 'p' if d == P else ('1' if d == 1 else '%d bits' % d.bit_length())
    say('  x_%-6d denom = %s   numerator %d bits' % (u, fac, abs(q.numerator).bit_length()))
    if d != 1:
        blocked.append((u, q))
say('\nblocking coordinates: %d' % len(blocked))

say('\n--- are the four conditions independent? residues of each numerator mod its denominator')
info = []
for u, q in blocked:
    d = q.denominator
    num = q.numerator
    res = num % d
    info.append((u, d, num, res))
    say('  x_%-6d  denom %d bits (== p: %s)  numerator mod denom = %d bits, zero: %s'
        % (u, d.bit_length(), d == P, res.bit_length(), res == 0))

say('\n--- relations among the residues (all taken mod p where the denominator is p)')
pres = [(u, r_) for u, d, nu, r_ in info if d == P]
for i in range(len(pres)):
    for j in range(i + 1, len(pres)):
        ui, ri = pres[i]
        uj, rj = pres[j]
        eq = (ri == rj)
        neg = ((ri + rj) % P == 0)
        try:
            ratio = (ri * pow(rj, -1, P)) % P
        except Exception:
            ratio = None
        say('  x_%d vs x_%d: equal=%s  negatives=%s  ratio mod p = %s'
            % (ui, uj, eq, neg, str(ratio)[:30] + '..' if ratio is not None else 'n/a'))

say('\n--- the non-p denominator')
for u, d, nu, r_ in info:
    if d != P:
        say('  x_%d denom = %s (%d bits)' % (u, d, d.bit_length()))
        say('     gcd(denom, p) = %d' % math.gcd(int(d), P))
        dd = int(d)
        small = []
        for q_ in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71,
                   73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 2264251, 11, 4787, 537773,
                   6672769, 5113045, 7376877, 14321763, 2169493, 4914685]:
            while dd % q_ == 0:
                small.append(q_)
                dd //= q_
        say('     small factors: %s  cofactor %d bits' % (sorted(set(small)), dd.bit_length()))
        say('     cofactor == p: %s   p | denom: %s' % (dd == P, int(d) % P == 0))

say('\n--- what the conditions mean structurally')
say('columns divisible by p:')
for u in Pv:
    cs = list(cols[u].values())
    say('  x_%-6d entries %s  all p-divisible: %s'
        % (u, [abs(c).bit_length() for c in cs], all(c % P == 0 for c in cs)))
say('DONE')
