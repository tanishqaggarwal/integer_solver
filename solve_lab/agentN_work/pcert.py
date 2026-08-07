"""A mod-p CERTIFICATE for the barrier, replacing a search with a one-line check.

On the exact zero-collateral variety (rank 14, complete 68-knob set) the region system is affine:
`M x = -b` over Z.  A necessary condition for integer solvability of a row subset is solvability of
the same subset MOD p.  If every 6-subset is inconsistent mod p, the "no six rows are integrally
zeroable" result is certified by a rank comparison rather than established by a search — and the
certificate is checkable in a line.

Both directions are reported: how many 6-subsets are inconsistent mod p, and (for honesty) how many
are consistent mod p yet still integrally unsolvable, i.e. where p is NOT the whole story.
"""
import os, sys, json, pickle, itertools

HERE = os.path.dirname(os.path.abspath(__file__))
sys.set_int_max_str_digits(2000000)
import zsolve

Pp = 115792089237316195423570985008687907853269984665640564039457584007908834671663
D = pickle.load(open(os.path.join(HERE, 'runs', 'solve68.pkl'), 'rb'))
M, b, Rl = D['M'], D['b'], D['Rl']
n = len(M[0])


def rank_mod(rows, ncol, p):
    A = [[x % p for x in r] for r in rows]
    r = 0
    for c in range(ncol):
        piv = None
        for i in range(r, len(A)):
            if A[i][c]:
                piv = i
                break
        if piv is None:
            continue
        A[r], A[piv] = A[piv], A[r]
        inv = pow(A[r][c], p - 2, p)
        A[r] = [(x * inv) % p for x in A[r]]
        for i in range(len(A)):
            if i != r and A[i][c]:
                f = A[i][c]
                A[i] = [(A[i][j] - f * A[r][j]) % p for j in range(ncol)]
        r += 1
    return r


def consistent_mod(idx, p):
    sub = [M[i] for i in idx]
    aug = [M[i] + [b[i]] for i in idx]
    return rank_mod(sub, n, p) == rank_mod(aug, n + 1, p)


print('region system %d x %d on the exact zero-collateral variety' % (len(M), n), flush=True)
for size in (5, 6, 7, 12):
    tot = incons = 0
    cons_examples = []
    for idx in itertools.combinations(range(12), size):
        tot += 1
        if not consistent_mod(idx, Pp):
            incons += 1
        elif len(cons_examples) < 5:
            cons_examples.append([Rl[i] for i in idx])
    print('subsets of size %-2d : %5d total, %5d INCONSISTENT mod p, %5d consistent mod p'
          % (size, tot, incons, tot - incons), flush=True)
    if cons_examples:
        print('     e.g. consistent mod p: %s' % cons_examples[:3], flush=True)

print('\nall-12 system: consistent mod p = %s' % consistent_mod(tuple(range(12)), Pp), flush=True)

# cross-check the certificate against the exact integer oracle on size 6
print('\ncross-check against the exact integer oracle (zsolve) on all 924 six-subsets:', flush=True)
agree = dis = 0
for idx in itertools.combinations(range(12), 6):
    cm = consistent_mod(idx, Pp)
    o, _, ex, _ = zsolve.max_zero_rows([M[i] for i in idx], [b[i] for i in idx], n, 6,
                                       node_cap=200000)
    intsolv = (o == 6)
    if (not cm) and (not intsolv):
        agree += 1
    elif cm and not intsolv:
        dis += 1
    elif cm and intsolv:
        agree += 1
    else:
        print('   IMPOSSIBLE: inconsistent mod p but integrally solvable, %s'
              % [Rl[i] for i in idx], flush=True)
print('   mod-p certificate settles %d of 924 ; %d need more than p' % (agree, dis), flush=True)

# is p the only prime that certifies?  try a few others for contrast
print('\ncontrast: the same certificate at other primes', flush=True)
for q in (2, 3, 5, 1000003, 2 ** 61 - 1):
    incons = sum(1 for idx in itertools.combinations(range(12), 6) if not consistent_mod(idx, q))
    print('   p = %-25s : %3d of 924 six-subsets inconsistent' % (q, incons), flush=True)

json.dump(dict(n=n), open(os.path.join(HERE, 'runs', 'pcert.json'), 'w'))
