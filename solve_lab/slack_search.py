#!/usr/bin/env python3
"""Search in SLACK-ACTIVE space. For each x_12779=1 activating bit, report:
  - v1 (plain forward, no slack) violated atoms  [cost of the bit itself]
  - v2 (slack-active, twist forced) violated atoms  [+ripple cost]
Then greedily add more bits to drive v2 violations down. The twist is always
satisfied by construction, so we minimize the SECONDARY atoms."""
import json, time
from confluent_eval5 import build5, make_forward
from slack_active import make_slack_solver, viol_atoms

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
    control = json.load(open('control_bits.json'))

    def slack_eval(bits):
        v1 = solve(list(bestval), bits)
        if v1[12779] != 1:
            return None, None, None
        frozen = {24026: v1[18274]-v1[35186], 27116: v1[17728]-v1[1642]}
        v2 = run(list(v1), frozen)
        return v1, v2, frozen

    # which single bits set x_12779=1?
    act = []
    for b in control:
        v = solve(list(bestval), [b])
        if v[12779] == 1: act.append(b)
    print(f"{len(act)} single bits set x_12779=1: {sorted(act)}", flush=True)

    # for each activator, v1 and v2 violation counts
    best_act = None
    for b in sorted(act):
        v1 = solve(list(bestval), [b])
        b1 = viol_atoms(A, v1)
        _, v2, fr = slack_eval([b])
        b2 = viol_atoms(A, v2)
        twist = (v2[9770]==v2[18274], v2[3183]==v2[17728])
        print(f"  bit {b}: v1_viol={len(b1)} {sorted(b1)[:6]}  v2_viol={len(b2)} twist={twist}", flush=True)
        if best_act is None or len(b2) < best_act[1]:
            best_act = (b, len(b2))
    print(f"\nbest single activator: bit {best_act[0]} with {best_act[1]} slack-active violations", flush=True)

    # greedy: from best activator, add bits to reduce v2 violations
    cur = [best_act[0]]
    _, v2, _ = slack_eval(cur)
    curv = len(viol_atoms(A, v2))
    print(f"\ngreedy from {cur}: {curv} violations", flush=True)
    improved = True; rounds = 0
    while improved and rounds < 12 and time.time()-t0 < 2400:
        improved = False; rounds += 1
        bestb = None; bestn = curv
        for b in control:
            if b in cur: continue
            r = slack_eval(cur+[b])
            if r[1] is None: continue
            n = len(viol_atoms(A, r[1]))
            if n < bestn: bestn = n; bestb = b
        if bestb is not None:
            cur.append(bestb); curv = bestn; improved = True
            print(f"  round {rounds}: +bit {bestb} -> {curv} violations ({time.time()-t0:.0f}s)", flush=True)
            if curv == 0:
                r = slack_eval(cur)
                json.dump({f"x_{i}": r[1][i] for i in range(len(r[1]))}, open('cand_slack_SOLVED.json','w'))
                print("  *** ALL ATOMS SATISFIED (slack-active) ***", flush=True); return
    print(f"\ngreedy plateau: {curv} violations with bits {sorted(cur)} ({time.time()-t0:.0f}s)", flush=True)

if __name__ == '__main__':
    main()
