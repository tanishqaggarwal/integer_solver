#!/usr/bin/env python3
"""Probe the multilinear degree of x_9770,x_3183 in the 22 bits (mod P, which is
artifact-free). Compute Mobius coefficients up to degree 3; report how many are
nonzero at each degree -> the polynomial's true degree. If low, we can tabulate
all 2^22 exactly and solve x_9770(A)=T1, x_3183(A)=T2."""
import json, time, itertools
from confluent_eval5 import build5, make_forward

BITS22 = [1782,1858,2795,2800,3483,5443,10652,19520,21188,21588,23634,26947,
          27512,29682,30104,30596,30658,30792,33251,37748,37885,38116]
P = (1 << 61) - 1

def main():
    t0 = time.time()
    A, kind, info, seq, bestval, ncyc = build5()
    order = json.load(open('eval_order.json'))['order']
    defset = set(v for v in kind if kind[v] != 'const')
    seq = [v for v in order if v in defset and v not in (9770, 3183)]
    seq += [v for v in (9770, 3183) if v in defset]
    seq += [v for v in defset if v not in set(order) and v not in (9770, 3183)]
    solve = make_forward(kind, info, seq, bestval, mod=P)
    bm = [x % P for x in bestval]

    # also Z solver for later exact checks
    def ev(S):
        val = solve(list(bm), list(S))
        return val[9770] % P, val[3183] % P

    f0 = ev([])
    print(f"f(0) = {f0}  ({time.time()-t0:.0f}s)", flush=True)
    fs = {b: ev([b]) for b in BITS22}
    # degree-1 coeffs
    c1 = {b: tuple((fs[b][i] - f0[i]) % P for i in range(2)) for b in BITS22}
    nz1 = [b for b in BITS22 if any(c1[b])]
    print(f"deg1 nonzero: {len(nz1)}/22 for the pair (9770,3183)", flush=True)

    # degree-2 coeffs
    c2 = {}
    for b, c in itertools.combinations(BITS22, 2):
        fbc = ev([b, c])
        val = tuple((fbc[i] - fs[b][i] - fs[c][i] + f0[i]) % P for i in range(2))
        c2[(b, c)] = val
    nz2_9770 = sum(1 for v in c2.values() if v[0])
    nz2_3183 = sum(1 for v in c2.values() if v[1])
    print(f"deg2 nonzero: 9770={nz2_9770}/231, 3183={nz2_3183}/231 ({time.time()-t0:.0f}s)", flush=True)

    # degree-3 coeffs (sample: all triples among first 10 bits, plus some spread)
    tri = list(itertools.combinations(BITS22[:10], 3))
    nz3_9770 = nz3_3183 = 0
    maxdeg3 = 0
    for b, c, d in tri:
        fbcd = ev([b, c, d])
        # mobius deg3 = f(bcd) - sum f(pairs) + sum f(singles) - f(0)
        val0 = (fbcd[0] - ev([b,c])[0] - ev([b,d])[0] - ev([c,d])[0]
                + fs[b][0] + fs[c][0] + fs[d][0] - f0[0]) % P
        val1 = (fbcd[1] - ev([b,c])[1] - ev([b,d])[1] - ev([c,d])[1]
                + fs[b][1] + fs[c][1] + fs[d][1] - f0[1]) % P
        if val0: nz3_9770 += 1
        if val1: nz3_3183 += 1
    print(f"deg3 (sample {len(tri)} triples): 9770 nz={nz3_9770}, 3183 nz={nz3_3183} ({time.time()-t0:.0f}s)", flush=True)

    # degree-4 sample among first 8
    quad = list(itertools.combinations(BITS22[:8], 4))
    nz4 = 0
    for Q in quad[:20]:
        # inclusion-exclusion for the 4-subset
        s = 0
        for r in range(len(Q)+1):
            for sub in itertools.combinations(Q, r):
                s = (s + (-1)**(len(Q)-r) * ev(list(sub))[0]) % P
        if s: nz4 += 1
    print(f"deg4 (sample {min(20,len(quad))}): 9770 nz={nz4} ({time.time()-t0:.0f}s)", flush=True)
    print("DONE", flush=True)

if __name__ == '__main__':
    main()
