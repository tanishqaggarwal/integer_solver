"""STEP 2 — are the blocking congruences a FUNCTION of the region's boundary data?

The region system is A z = b.  A is fixed by the instance (its entries are p, ±1, and the
literal multipliers).  b is built from the region's boundary data: the residual constants
const_a, which are affine in the externally-determined quantities
    J  = x_7075*x_8731,  K1 = x_7068 - x_2099,  K2 = 5113045*x_7075*x_9118,
    L  = x_4432 - x_19964   (the a23618 constant)
A configuration scan varies exactly those.  So: perturb each const_a by +1 and measure the
induced change in each blocking residue mod its modulus.  If the map is affine and its image
is all of (Z/p)^4 x (Z/2458959)^2, the conditions can be INVERTED — the scan becomes a
targeted construction, not a sample.  If the image is a proper subgroup, the scan is a
measurement and the rate is the subgroup index.
"""
import sys, json, math
from fractions import Fraction
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentO_work')
import regiongrow as G, harness as H

OD = '/home/user/integer_solver/solve_lab/agentO_work'
LOG = open(OD + '/runs/structure.log', 'w', buffering=1)
P = 115792089237316195423570985008687907853269984665640564039457584007908834671663
SMALL = 2458959


def say(*a):
    print(*a, file=LOG)


R = G.R0 + [23618]
Pv = G.private_vars(R)
const0, cols = G.build_model(R, Pv, G.V0)
Eqs, rows0 = G.eq_system(R, Pv, const0, cols)
A = [[Fraction(rows0[e][0].get(u, 0)) for u in Pv] for e in Eqs]


def solve_rat(const):
    _, rows = G.eq_system(R, Pv, const, cols)
    b = [Fraction(rows[e][1]) for e in Eqs]
    M = [A[i][:] + [b[i]] for i in range(len(A))]
    n, m = len(M), len(Pv)
    piv, r = [], 0
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
    if any(M[i][m] != 0 for i in range(r, n)):
        return None
    return {Pv[c]: M[i][m] for i, c in enumerate(piv)}


BLOCK = [(1329, P), (9413, P), (10903, P), (17325, P), (642, SMALL), (17325, SMALL)]


def residues(z):
    out = []
    for u, mod in BLOCK:
        q = z[u]
        # q = num/den; the obstruction is  num * den^-1  mod  mod   (den coprime to mod part)
        num, den = q.numerator, q.denominator
        d = den
        # strip the part of the denominator equal to mod
        k = 0
        while d % mod == 0:
            d //= mod
            k += 1
        if k == 0:
            out.append(0)
            continue
        try:
            inv = pow(d % mod, -1, mod)
        except ValueError:
            out.append(None)
            continue
        out.append((num * inv) % mod)
    return out


z0 = solve_rat(const0)
r0 = residues(z0)
say('blocking conditions (coordinate, modulus) and their residues at the witness:')
for (u, mod), rv in zip(BLOCK, r0):
    say('   x_%-6d mod %s : residue %s  (zero: %s)'
        % (u, 'p' if mod == P else str(mod), (str(rv)[:28] + '..') if rv and rv > 10 ** 12 else rv, rv == 0))

say('\n--- sensitivity: perturb each region boundary constant by +1')
rowsout = []
for a in sorted(const0):
    c = dict(const0)
    c[a] = const0[a] + 1
    z = solve_rat(c)
    if z is None:
        say('   const a%-6d: perturbation makes the system rationally inconsistent' % a)
        continue
    r1 = residues(z)
    dv = []
    for (u, mod), x0, x1 in zip(BLOCK, r0, r1):
        dv.append(None if (x0 is None or x1 is None) else (x1 - x0) % mod)
    rowsout.append((a, dv))
    say('   const a%-6d -> d(residue) = %s'
        % (a, [('0' if d == 0 else (str(d)[:14] + '..' if d and d > 10 ** 12 else str(d))) for d in dv]))

say('\n--- linearity check: perturb by +2, expect exactly twice the +1 delta')
lin = True
for a in sorted(const0):
    c = dict(const0)
    c[a] = const0[a] + 2
    z = solve_rat(c)
    if z is None:
        continue
    r2 = residues(z)
    d1 = dict(rowsout).get(a)
    if d1 is None:
        continue
    for i, ((u, mod), x0, x2) in enumerate(zip(BLOCK, r0, r2)):
        if x0 is None or x2 is None:
            continue
        if (x2 - x0) % mod != (2 * d1[i]) % mod:
            lin = False
            say('   a%d coordinate %d NOT linear' % (a, i))
say('   all perturbations linear: %s' % lin)

say('\n--- rank of the sensitivity map mod p (4 p-conditions) and mod %d (2 conditions)' % SMALL)


def rank_mod(rowsM, mod):
    Mx = [row[:] for row in rowsM]
    n = len(Mx)
    m = len(Mx[0]) if n else 0
    rk = 0
    for c in range(m):
        p_ = None
        for i in range(rk, n):
            if Mx[i][c] % mod:
                try:
                    pow(Mx[i][c] % mod, -1, mod)
                except ValueError:
                    continue
                p_ = i
                break
        if p_ is None:
            continue
        Mx[rk], Mx[p_] = Mx[p_], Mx[rk]
        inv = pow(Mx[rk][c] % mod, -1, mod)
        Mx[rk] = [(x * inv) % mod for x in Mx[rk]]
        for i in range(n):
            if i != rk and Mx[i][c] % mod:
                f = Mx[i][c] % mod
                Mx[i] = [(Mx[i][k] - f * Mx[rk][k]) % mod for k in range(m)]
        rk += 1
    return rk


Mp = [[rowsout[j][1][i] or 0 for i in range(4)] for j in range(len(rowsout))]
Ms = [[rowsout[j][1][i] or 0 for i in range(4, 6)] for j in range(len(rowsout))]
say('   boundary constants available: %d' % len(rowsout))
say('   rank of sensitivity mod p     = %d  (need 4 to invert all p-conditions)' % rank_mod(Mp, P))
say('   rank of sensitivity mod %d = %d  (need 2)' % (SMALL, rank_mod(Ms, SMALL)))
json.dump({'residues': [str(x) for x in r0],
           'sens': {str(a): [str(x) for x in d] for a, d in rowsout}},
          open(OD + '/structure.json', 'w'))
say('DONE')
