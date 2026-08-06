#!/usr/bin/env python3
"""Test linearity (mod P) of the HUGE residue cores x_26517, x_22190 in the 233
bits, and characterize the small corrections x_34150, x_24424. Trace their gates.
If the huge cores are linear subset-sums, the target eqs become lattice-solvable."""
import json
from confluent_eval5 import build5, make_forward
from propagate import atom_vars

BITS22 = set([1782,1858,2795,2800,3483,5443,10652,19520,21188,21588,23634,26947,
          27512,29682,30104,30596,30658,30792,33251,37748,37885,38116])
P = (1 << 61) - 1

def main():
    A, kind, info, seq0, bestval, ncyc = build5()
    order = json.load(open('eval_order.json'))['order']
    defset = set(v for v in kind if kind[v] != 'const')
    seq = [v for v in order if v in defset and v not in (9770, 3183)]
    seq += [v for v in (9770, 3183) if v in defset]
    seq += [v for v in defset if v not in set(order) and v not in (9770, 3183)]
    solve = make_forward(kind, info, seq, bestval)
    solveP = make_forward(kind, info, seq, bestval, mod=P)
    control = json.load(open('control_bits.json'))
    bits233 = [b for b in control if b not in BITS22]
    base = solve(list(bestval), [])
    bm = [x % P for x in bestval]

    def cls(v):
        x = abs(base[v]);
        return '0' if x==0 else (f'small({base[v]})' if x<10**8 else f'HUGE({len(str(x))}d)')
    for w in (26517, 34150, 22190, 24424):
        k = kind[w]
        s = f"kind={k}"
        if k == 'gate':
            coef, terms = info[w]
            s += f" coef={coef} nterms={len(terms)}: " + " ".join(
                ('*'.join(f'x_{x}[{cls(x)}]' for x in m)) for c,m in terms[:4])
        elif k == 'div':
            c,u,rest = info[w]; s += f" div u=x_{u}[{cls(u)}] rest={len(rest)}"
        elif k == 'load':
            s += " LOAD"
        print(f"x_{w} [{cls(w)}] {s}")

    # linearity of x_26517, x_22190 over 233 bits, mod P
    b0 = solveP(list(bm), [])
    WCHK = [26517, 22190, 34150, 24424]
    single = {b: solveP(list(bm), [b]) for b in bits233}
    def lintest(w, sizes, seedmul):
        out = []
        for si, sz in enumerate(sizes):
            S = sorted(set(bits233[(i*seedmul+si*131+3) % len(bits233)] for i in range(sz)))
            v = solveP(list(bm), S)
            actual = v[w]
            pred = (b0[w] + sum(single[b][w]-b0[w] for b in S)) % P
            out.append('LIN' if actual==pred else 'NL')
        print(f"  x_{w}: {out} at {sizes}", flush=True)
    print("linearity over 233 bits (mod P):")
    for w in WCHK: lintest(w, [3,8,30,80,150,233], 91)

    # how many bits move x_34150 / x_24424, and their value ranges (Z, sparse)
    for w in (34150, 24424):
        movers = [b for b in bits233 if solve(list(bestval),[b])[w] != base[w]]
        print(f"  x_{w}: moved by {len(movers)} of 233 bits; base={base[w]}")

if __name__ == '__main__':
    main()
