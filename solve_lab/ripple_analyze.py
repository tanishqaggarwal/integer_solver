#!/usr/bin/env python3
"""Analyze the atoms broken by slack activation (bit=1858 gave 18). For each, show
the polynomial, which vars are 'perturbed' (differ from forward-eval baseline), and
whether the atom has a slack var we could still choose. Goal: understand the ripple
and find a consistent slack-active assignment."""
import json, time
from confluent_eval5 import build5, make_forward
from slack_active import make_slack_solver, viol_atoms
from propagate import atom_vars

def fmt(poly, control):
    parts = []
    for m, c in sorted(poly.items(), key=lambda kv: (-len(kv[0]), kv[0])):
        cc = str(c) if abs(c) < 10**9 else f'H{len(str(abs(c)))}'
        mon = '*'.join('x'+str(x)+('#' if x in control else '') for x in m) if m else '1'
        parts.append(f'{cc}*{mon}')
    return ' + '.join(parts)

def main():
    t0 = time.time()
    A, kind, info, seq0, bestval, ncyc = build5()
    order = json.load(open('eval_order.json'))['order']
    defset = set(v for v in kind if kind[v] != 'const')
    seq = [v for v in order if v in defset and v not in (9770, 3183)]
    seq += [v for v in (9770, 3183) if v in defset]
    seq += [v for v in defset if v not in set(order) and v not in (9770, 3183)]
    solve = make_forward(kind, info, seq, bestval)
    run, seq2 = make_slack_solver(kind, info, seq, bestval)
    control = set(json.load(open('control_bits.json')))

    bit = 1858
    v1 = solve(list(bestval), [bit])
    frozen = {24026: v1[18274]-v1[35186], 27116: v1[17728]-v1[1642]}
    v2 = run(list(v1), frozen)
    bad = viol_atoms(A, v2)
    print(f"bit={bit}: {len(bad)} broken atoms; x_24026={frozen[24026]}, x_27116={frozen[27116]}", flush=True)
    # perturbed vars: differ between v1 (baseline forward) and v2 (slack-active)
    pert = set(v for v in range(len(v1)) if v1[v] != v2[v])
    print(f"perturbed vars (v1!=v2): {len(pert)}", flush=True)
    for a in sorted(bad):
        poly = A[a]
        s = 0
        for m, c in poly.items():
            t = c
            for x in m: t *= v2[x]
            s += t
        pv = sorted(atom_vars(poly) & pert)
        # find vars in atom that are NOT perturbed and NOT control (could be free knobs)
        av = atom_vars(poly)
        print(f"\n  a{a} [nv={len(av)} resid={s}]  perturbed-in-atom={pv}")
        print(f"     {fmt(poly, control)[:150]}")
    print(f"\ndone ({time.time()-t0:.0f}s)", flush=True)

if __name__ == '__main__':
    main()
