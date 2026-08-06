#!/usr/bin/env python3
"""Extract exact structure of the twist/cascade checks 1817,44271,30378,40782."""
import json
from confluent_eval5 import build5, make_forward
from propagate import NVARS, atom_vars

def main():
    A, kind, info, seq, bestval, ncyc = build5()
    order = json.load(open('eval_order.json'))['order']
    defset = set(v for v in kind if kind[v] != 'const')
    seq = [v for v in order if v in defset and v not in (9770, 3183)]
    seq += [v for v in (9770, 3183) if v in defset]
    seq += [v for v in defset if v not in set(order) and v not in (9770, 3183)]
    solve = make_forward(kind, info, seq, bestval)
    base = solve(list(bestval), [])
    for a in (1817, 44271, 30378):
        print(f"\n=== atom {a} (resid0={sum_resid(A[a], base)}) ===")
        for m, c in sorted(A[a].items(), key=lambda kv: kv[0]):
            vv = ' * '.join(f'x_{x}' for x in m) if m else '1'
            print(f"   {c:+d} * {vv}   [x_{m[0]}={base[m[0]]}]" if len(m)==1 else f"   {c:+d} * {vv}")
    # atom 40782: which watch-vars / control appear
    va = atom_vars(A[40782])
    control = set(json.load(open('control_bits.json')))
    BITS22 = {1782,1858,2795,2800,3483,5443,10652,19520,21188,21588,23634,26947,27512,29682,30104,30596,30658,30792,33251,37748,37885,38116}
    print(f"\n=== atom 40782: {len(va)} vars ===")
    print("  contains x_9770?", 9770 in va, " x_3183?", 3183 in va, " x_18274?", 18274 in va, " x_17728?", 17728 in va)
    print("  control bits in atom:", sorted(va & control))
    print("  vars:", sorted(va))
    # kind of the watch vars
    for w in (9770, 3183, 18274, 17728):
        print(f"  kind[x_{w}] = {kind.get(w)}")

def sum_resid(poly, val):
    s = 0
    for m, c in poly.items():
        t = c
        for x in m: t *= val[x]
        s += t
    return s

if __name__ == '__main__':
    main()
