#!/usr/bin/env python3
"""a40782=0 is implied by the FULL twist (x_9770=x_18274 & x_17728=x_3183) per
test_40782. My slack-active state satisfies both twist halves yet a40782 is broken.
Diagnose: evaluate Q40782 (deg-2 root) at the slack-active state, list which terms
are nonzero, and check the auxiliary twist vars (x_9982=x_9897*x_12518, x_26977=
x_20510*x_31302) that a40782 also needs at 0. Do the same for a39550. This tells us
EXACTLY the residual constraints to satisfy."""
import json, time
from confluent_eval5 import build5, make_forward
from slack_active import make_slack_solver
from check_square import try_sqrt
from propagate import atom_vars

def main():
    t0 = time.time()
    A0, kind, info, seq0, bestval, ncyc = build5()
    order = json.load(open('eval_order.json'))['order']
    defset = set(v for v in kind if kind[v] != 'const')
    seq = [v for v in order if v in defset and v not in (9770, 3183)]
    seq += [v for v in (9770, 3183) if v in defset]
    seq += [v for v in defset if v not in set(order) and v not in (9770, 3183)]
    solve = make_forward(kind, info, seq, bestval)
    run, seq2 = make_slack_solver(kind, info, seq, bestval)

    v1 = solve(list(bestval), [1858])
    frozen = {24026: v1[18274]-v1[35186], 27116: v1[17728]-v1[1642]}
    val = run(list(v1), frozen)
    print(f"twist: x_9770==x_18274 {val[9770]==val[18274]}, x_3183==x_17728 {val[3183]==val[17728]}", flush=True)
    aux = {9982:'x_9897*x_12518', 26977:'x_20510*x_31302', 9897:'', 12518:'', 20510:'', 31302:''}
    for v in aux:
        print(f"  x_{v} = {val[v]}   {aux[v]}", flush=True)

    for a in (40782, 39550):
        Q = try_sqrt(A0[a])
        # evaluate Q, list nonzero terms grouped by variable
        qval = 0; nzterms = []
        for m, c in Q.items():
            t = c
            for x in m: t *= val[x]
            if t != 0: nzterms.append((m, t))
            qval += t
        print(f"\na{a}: Q residual = {qval}", flush=True)
        print(f"  Q has {len(Q)} terms, {len(nzterms)} nonzero at slack-active state", flush=True)
        # which variables appear in the nonzero terms?
        nzvars = set()
        for m, t in nzterms: nzvars.update(m)
        print(f"  vars in nonzero terms: {sorted(nzvars)}", flush=True)
        # value of each such var
        for v in sorted(nzvars):
            print(f"    x_{v} = {val[v]}", flush=True)
    print(f"\ndone ({time.time()-t0:.0f}s)", flush=True)

if __name__ == '__main__':
    main()
