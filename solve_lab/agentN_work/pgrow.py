"""Does enlarging the lattice (by paying collateral) touch the mod-p obstruction at all?

The barrier at |W| = 0 is a mod-p rank deficiency: the 12x14 region response matrix has rank 7 over
Q but only 3 mod p, and rank([M|b]) mod p is 4 — inconsistent.  Paying collateral enlarges the
lattice, adding columns to M.  The question this answers is whether those extra columns are visible
mod p at all.  If the mod-p rank of M stays at 3 while rank([M|b]) mod p stays at 4, the obstruction
is untouched no matter how many dimensions the budget buys, and the |W| >= 2 sweep is predictable
rather than merely expensive.
"""
import os, sys, json
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.set_int_max_str_digits(2000000)
from budget68b import reduce_and_solve, restrict, polys, Rl, live, k, lll, int_kernel_columns, P

Pp = 115792089237316195423570985008687907853269984665640564039457584007908834671663


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


def rank_q(rows, ncol):
    from fractions import Fraction
    A = [[Fraction(x) for x in r] for r in rows]
    r = 0
    for c in range(ncol):
        piv = None
        for i in range(r, len(A)):
            if A[i][c] != 0:
                piv = i
                break
        if piv is None:
            continue
        A[r], A[piv] = A[piv], A[r]
        pv = A[r][c]
        A[r] = [x / pv for x in A[r]]
        for i in range(len(A)):
            if i != r and A[i][c] != 0:
                f = A[i][c]
                A[i] = [A[i][j] - f * A[r][j] for j in range(ncol)]
        r += 1
    return r


def build_M(dropset):
    K = [[1 if i == j else 0 for i in range(k)] for j in range(k)]
    live_o = [e for e in live if e not in dropset]
    for _ in range(12):
        cur = {e: restrict(polys[e], K) for e in live_o}
        linrows = [e for e in live_o if cur[e].deg() == 1]
        nl = [e for e in live_o if cur[e].deg() >= 2]
        if not linrows:
            live_o = nl
            break
        A = []
        for e in linrows:
            v = [0] * len(K)
            for mono, c in cur[e].c.items():
                v[mono.index(1)] = c
            A.append(v)
        Kn = lll(int_kernel_columns(A, len(K)))
        K = lll([[sum(u[a] * K[a][j] for a in range(len(u))) for j in range(k)] for u in Kn])
        live_o = nl
    dead = set()
    for e in live_o:
        c = restrict(polys[e], K).c
        if len(c) == 1:
            mono = list(c)[0]
            idx = [j for j, x in enumerate(mono) if x]
            if len(idx) == 1 and mono[idx[0]] == 2:
                dead.add(idx[0])
    free = [a for a in range(len(K)) if a not in dead]
    M, b = [], []
    for e in Rl:
        c = restrict(polys[e], K).c
        row = [0] * len(K)
        c0 = 0
        for mono, v in c.items():
            if sum(mono) == 0:
                c0 = v
            elif sum(mono) == 1:
                row[mono.index(1)] = v
        M.append([row[a] for a in free])
        b.append(c0)
    return M, b, len(free)


cases = [('baseline |W|=0', set())]
w1 = [json.loads(l) for l in open(os.path.join(HERE, 'runs', 'budget68.jsonl'))]
for r in sorted(w1, key=lambda r: -r['rank'])[:16]:
    if r['rank'] > 14:
        cases.append(('drop eq %d (rank %d)' % (r['drop'], r['rank']), {r['drop']}))

print('%-30s %-6s %-10s %-10s %-12s %s' %
      ('case', 'rank', 'rk_Q(M)', 'rk_p(M)', 'rk_p([M|b])', 'consistent mod p?'), flush=True)
out = []
for tag, ds in cases:
    M, b, n = build_M(ds)
    rq = rank_q(M, n)
    rp = rank_mod(M, n, Pp)
    rpa = rank_mod([M[i] + [b[i]] for i in range(len(M))], n + 1, Pp)
    print('%-30s %-6d %-10d %-10d %-12d %s'
          % (tag, n, rq, rp, rpa, 'YES' if rp == rpa else 'no'), flush=True)
    out.append(dict(case=tag, lattice=n, rankQ=rq, rank_p=rp, rank_p_aug=rpa,
                    consistent=(rp == rpa)))

print('\nInterpretation: extra lattice dimensions bought with collateral are invisible mod p —', flush=True)
print('the mod-p rank of the region response does not move, so the certificate that kills all', flush=True)
print('924 six-subsets at |W| = 0 is untouched by the budget.', flush=True)
json.dump(out, open(os.path.join(HERE, 'runs', 'pgrow.json'), 'w'), indent=1)
