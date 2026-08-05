#!/usr/bin/env python3
"""Probe multilinear degree (mod P) of x_8821 (shared denom), x_6773, x_17233
(numerators) in the 233 bits. Low degree => algebraic attack on
NUM - N*DEN = 0 becomes tractable. Also trace their defining gates."""
import json, time, itertools
from confluent_eval5 import build5, make_forward

BITS22 = set([1782,1858,2795,2800,3483,5443,10652,19520,21188,21588,23634,26947,
          27512,29682,30104,30596,30658,30792,33251,37748,37885,38116])
P = (1 << 61) - 1
WIRES = [8821, 6773, 17233, 18274]

def main():
    t0 = time.time()
    A, kind, info, seq0, bestval, ncyc = build5()
    order = json.load(open('eval_order.json'))['order']
    defset = set(v for v in kind if kind[v] != 'const')
    seq = [v for v in order if v in defset and v not in (9770, 3183)]
    seq += [v for v in (9770, 3183) if v in defset]
    seq += [v for v in defset if v not in set(order) and v not in (9770, 3183)]
    control = json.load(open('control_bits.json'))
    bits233 = [b for b in control if b not in BITS22]
    solveP = make_forward(kind, info, seq, bestval, mod=P)
    bm = [x % P for x in bestval]

    for w in WIRES:
        print(f"kind[x_{w}]={kind[w]}", info.get(w) if kind[w] != 'gate' else f"gate coef={info[w][0]} nterms={len(info[w][1])}")

    def ev(S):
        val = solveP(list(bm), list(S))
        return {w: val[w] for w in WIRES}

    f0 = ev([])
    # use a fixed small set of bits for degree probing
    probe = bits233[:14]
    fs = {b: ev([b]) for b in probe}
    # deg2 over probe pairs
    c2nz = {w: 0 for w in WIRES}
    pairvals = {}
    for b, c in itertools.combinations(probe, 2):
        v = ev([b, c]); pairvals[(b, c)] = v
        for w in WIRES:
            m2 = (v[w] - fs[b][w] - fs[c][w] + f0[w]) % P
            if m2: c2nz[w] += 1
    npairs = len(list(itertools.combinations(probe, 2)))
    print(f"deg2 nonzero / {npairs}: " + ", ".join(f"x_{w}={c2nz[w]}" for w in WIRES), flush=True)
    # deg3 over probe triples
    c3nz = {w: 0 for w in WIRES}
    ntri = 0
    for b, c, d in itertools.combinations(probe, 3):
        ntri += 1
        vv = ev([b, c, d])
        for w in WIRES:
            m3 = (vv[w] - pairvals[(b,c)][w] - pairvals[(b,d)][w] - pairvals[(c,d)][w]
                  + fs[b][w] + fs[c][w] + fs[d][w] - f0[w]) % P
            if m3: c3nz[w] += 1
    print(f"deg3 nonzero / {ntri}: " + ", ".join(f"x_{w}={c3nz[w]}" for w in WIRES), flush=True)
    # deg4 over probe quads (first 10)
    c4nz = {w: 0 for w in WIRES}
    nq = 0
    for Q in itertools.combinations(probe[:9], 4):
        nq += 1
        s = {w: 0 for w in WIRES}
        for r in range(5):
            for sub in itertools.combinations(Q, r):
                vv = ev(list(sub))
                for w in WIRES:
                    s[w] = (s[w] + (-1)**(4-r) * vv[w]) % P
        for w in WIRES:
            if s[w]: c4nz[w] += 1
    print(f"deg4 nonzero / {nq}: " + ", ".join(f"x_{w}={c4nz[w]}" for w in WIRES), flush=True)
    print(f"done ({time.time()-t0:.0f}s)", flush=True)

if __name__ == '__main__':
    main()
