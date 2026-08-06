#!/usr/bin/env python3
"""x_18274 = NUM/DEN (division wire). If NUM(B), DEN(B) are LINEAR (subset-sum)
in the 233 bits, then for fixed target N1 the equation NUM - N1*DEN = 0 is linear
in B -> lattice-solvable. Test linearity of NUM,DEN (and the analogous wires for
x_17728) in the 233 bits, mod P (artifact-free)."""
import json, time
from confluent_eval5 import build5, make_forward

BITS22 = set([1782,1858,2795,2800,3483,5443,10652,19520,21188,21588,23634,26947,
          27512,29682,30104,30596,30658,30792,33251,37748,37885,38116])
P = (1 << 61) - 1

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

    for w in (18274, 17728):
        print(f"kind[x_{w}]={kind[w]}, info={info[w] if kind[w]!='gate' else ('gate',)}")
    # x_18274 div: numerator = -sum(rest), denom = c*x_u
    c18, u18, rest18 = info[18274]
    c17, u17, rest17 = info[17728]
    print(f"x_18274: DEN = {c18} * x_{u18}; NUM = -(sum of {len(rest18)} rest terms)")
    print(f"x_17728: DEN = {c17} * x_{u17}; NUM = -(sum of {len(rest17)} rest terms)")

    solveP = make_forward(kind, info, seq, bestval, mod=P)
    bm = [x % P for x in bestval]

    def numden(val, c, u, rest):
        num = 0
        for cc, m in rest:
            t = cc % P
            for x in m: t = (t*val[x]) % P
            num = (num + t) % P
        num = (-num) % P
        den = (c * val[u]) % P
        return num, den

    def eval_wires(S):
        val = solveP(list(bm), list(S))
        n18, d18 = numden(val, c18, u18, rest18)
        n17, d17 = numden(val, c17, u17, rest17)
        return {'N18': n18, 'D18': d18, 'N17': n17, 'D17': d17,
                'x18274': val[18274], 'x17728': val[17728], 'u18': val[u18], 'u17': val[u17]}

    base = eval_wires([])
    single = {b: eval_wires([b]) for b in bits233}
    print(f"got singles ({time.time()-t0:.0f}s)", flush=True)

    def lintest(key, sizes, seedmul):
        results = []
        for si, sz in enumerate(sizes):
            idx = [((i*seedmul + si*7919 + 3) % len(bits233)) for i in range(sz)]
            S = sorted(set(bits233[i] for i in idx))
            actual = eval_wires(S)[key]
            pred = (base[key] + sum((single[b][key] - base[key]) for b in S)) % P
            ok = (actual == pred)
            results.append(ok)
        return results

    sizes = [3, 8, 30, 80, 150, 233]
    for key in ('N18', 'D18', 'N17', 'D17', 'u18', 'u17', 'x18274', 'x17728'):
        r = lintest(key, sizes, 101)
        print(f"  {key}: linear? {['LIN' if x else 'NL' for x in r]} at sizes {sizes}", flush=True)
    print(f"done ({time.time()-t0:.0f}s)", flush=True)

if __name__ == '__main__':
    main()
