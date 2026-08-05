#!/usr/bin/env python3
"""With x_12779=1 (x_14402=0), a1813 (x_14402*x_24026=321447*x_38215=0) is 0=0, freeing x_24026;
a1660 makes x_3368=x_24026, so x_9770=x_35186+x_24026 can hit any target. Same for x_27116/x_3183.
CHECK: are x_24026, x_27116 constrained by any OTHER atoms? If they only appear in a1660/a1813
(and a1815/a27976), the twist is FREELY satisfiable once x_12779=1. Also find which bits set
x_12779=1."""
import json
from propagate import load_atoms, atom_vars
from confluent_eval5 import build5, make_forward

def main():
    atoms = load_atoms()
    best = json.load(open('best/best_partial_39019.json')); bv = {int(k[2:]): v for k, v in best.items()}
    control = set(json.load(open('control_bits.json')))
    def val(v): return bv.get(v, 0)
    # all atoms containing x_24026 / x_27116 / x_12779 / x_14402
    for target in (24026, 27116, 12779, 14402, 38215, 29437):
        ats = [a for a in range(len(atoms)) if target in atom_vars(atoms[a])]
        print(f"\n=== x_{target} appears in {len(ats)} atoms ===")
        for a in ats:
            poly = atoms[a]
            tt = []
            for m, c in sorted(poly.items(), key=lambda kv: (-len(kv[0]), kv[0])):
                cc = str(c) if abs(c) < 10**8 else f'H{len(str(abs(c)))}'
                tt.append(f'{cc}*' + ('*'.join('x'+str(x)+('#' if x in control else '') for x in m) if m else '1'))
            # residual at best
            s = 0
            for m, c in poly.items():
                t = c
                for x in m: t *= val(x)
                s += t
            print(f"  a{a} (resid@best={s}): " + ' + '.join(tt)[:130])

    # find bits that make x_12779 = 1
    A, kind, info, seq0, bestval, ncyc = build5()
    order = json.load(open('eval_order.json'))['order']
    defset = set(v for v in kind if kind[v] != 'const')
    seq = [v for v in order if v in defset and v not in (9770, 3183)]
    seq += [v for v in (9770, 3183) if v in defset]
    seq += [v for v in defset if v not in set(order) and v not in (9770, 3183)]
    solve = make_forward(kind, info, seq, bestval)
    base = solve(list(bestval), [])
    print(f"\nx_12779 at all-0: {base[12779]}; x_14402: {base[14402]}")
    got1 = []
    for b in control:
        v = solve(list(bestval), [b])[12779]
        if v == 1: got1.append(b)
    print(f"single bits giving x_12779==1: {got1[:20]}")

if __name__ == '__main__':
    main()
