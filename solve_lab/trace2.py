#!/usr/bin/env python3
"""Trace x_6773 and x_17233 (the two numerators). Classify each wire in their
immediate cones as HUGE (residue) vs SMALL (selector/count). Test whether each
numerator is a SUBSET-SUM of huge residues gated by the 233 bits (i.e. linear in
the bits when the SELECTOR/denominator wires are held fixed)."""
import json, time
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
    control = json.load(open('control_bits.json'))
    bits233 = [b for b in control if b not in BITS22]
    base = solve(list(bestval), [])

    def cls(v):
        x = abs(base[v])
        if x == 0: return '0'
        if x < 10**8: return f'small({base[v]})'
        return f'HUGE({len(str(x))}d)'

    for w in (6773, 17233, 8821):
        print(f"\n=== x_{w} (kind={kind[w]}, val={cls(w)}) ===")
        if kind[w] == 'gate':
            coef, terms = info[w]
            print(f"  x_{w} = -(1/{coef}) * sum of {len(terms)} terms:")
            for c, m in terms:
                desc = ' * '.join(f'x_{x}[{cls(x)}]' for x in m)
                print(f"     {c:+d} * {desc}")
        elif kind[w] == 'div':
            c, u, rest = info[w]
            print(f"  x_{w} = -(sum rest)/({c}*x_{u}[{cls(u)}]); rest={len(rest)} terms")
            for cc, m in rest[:6]:
                desc = ' * '.join(f'x_{x}[{cls(x)}]' for x in m)
                print(f"     {cc:+d} * {desc}")

    # test: is x_6773 linear in the 233 bits when the 18 x_8821-bits are held at 0?
    rp = json.load(open('residue_pool.json'))
    x8821_bits = [int(b) for b in rp['x8821_weights']]
    print(f"\nx_8821 depends on bits: {sorted(x8821_bits)}")
    freebits = [b for b in bits233 if b not in set(x8821_bits)]
    print(f"233-side bits NOT feeding x_8821: {len(freebits)}")

    # linearity of x_6773 over freebits (x_8821 stays 1), mod P
    solveP = make_forward(kind, info, seq, bestval, mod=P)
    bm = [x % P for x in bestval]
    b0 = solveP(list(bm), [])
    single = {}
    for b in freebits:
        v = solveP(list(bm), [b]); single[b] = (v[6773], v[17233], v[8821])
    def lintest(idx, name, sizes, seedmul):
        for si, sz in enumerate(sizes):
            S = sorted(set(freebits[(i*seedmul+si*131+3) % len(freebits)] for i in range(sz)))
            v = solveP(list(bm), S)
            actual = v[[6773,17233,8821][idx]]
            pred = (b0[[6773,17233,8821][idx]] + sum(single[b][idx]-b0[[6773,17233,8821][idx]] for b in S)) % P
            print(f"    {name} |S|={sz} (x_8821={v[8821]}): {'LIN' if actual==pred else 'NL'}", flush=True)
    print("with x_8821-feeding bits held 0:")
    lintest(0, 'x_6773', [3,8,30,80,150,len([1 for _ in range(0)])+200 if False else 200], 71)
    lintest(1, 'x_17233', [3,8,30,80,150,200], 73)

if __name__ == '__main__':
    main()
