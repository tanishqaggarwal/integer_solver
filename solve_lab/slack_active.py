#!/usr/bin/env python3
"""SLACK-ACTIVE evaluator. Forward-eval forces x_24026=x_27116=0 (div wires a1813/
a1815 divide by x_14402=1-x_12779, =0 when x_12779=1). This makes the twist look
rigid. Here we RE-ORIENT: freeze x_24026,x_27116 as exogenous inputs, set
  x_24026 := x_18274 - x_35186   (=> x_9770 = x_35186 + x_12779*x_24026 = x_18274)
  x_27116 := x_17728 - x_1642     (=> x_3183 = x_1642 + x_12779*x_27116 = x_17728)
after computing the side-values in pass 1, with x_12779=1. Then re-run forward-eval
with those frozen, and COUNT violated atoms. This reaches the slack-active witness
state that plain forward-eval cannot represent."""
import json, time, sys
from confluent_eval5 import build5
from propagate import atom_vars, NVARS

def viol_atoms(A, val):
    bad = []
    for a, poly in enumerate(A):
        s = 0
        for m, c in poly.items():
            t = c
            for x in m: t *= val[x]
            s += t
        if s: bad.append(a)
    return bad

def make_slack_solver(kind, info, seq, bestval):
    FREEZE = {24026, 27116}
    seq2 = [v for v in seq if v not in FREEZE]
    def run(val, frozen):
        for v in FREEZE:
            val[v] = frozen.get(v, val[v])
        for v in seq2:
            k = kind[v]
            if k == 'gate':
                coef, terms = info[v]; rs = 0
                for c, m in terms:
                    t = c
                    for x in m: t *= val[x]
                    rs += t
                if coef and (-rs) % coef == 0: val[v] = (-rs)//coef
            elif k == 'load':
                bit, cbx, lt = info[v]
                if val[bit] == 0: val[v] = 0
                else:
                    rest = 0
                    for c, m in lt:
                        t = c
                        for x in m: t *= (1 if x == bit else val[x])
                        rest += t
                    num = -rest; den = cbx*val[bit]
                    if den and num % den == 0: val[v] = num//den
            elif k == 'div':
                c, u, rest = info[v]; rs = 0
                for cc, m in rest:
                    t = cc
                    for x in m: t *= val[x]
                    rs += t
                den = c*val[u]
                if den and (-rs) % den == 0: val[v] = (-rs)//den
                elif den == 0: val[v] = 0
        return val
    return run, seq2

def main():
    t0 = time.time()
    A, kind, info, seq0, bestval, ncyc = build5()
    order = json.load(open('eval_order.json'))['order']
    defset = set(v for v in kind if kind[v] != 'const')
    seq = [v for v in order if v in defset and v not in (9770, 3183)]
    seq += [v for v in (9770, 3183) if v in defset]
    seq += [v for v in defset if v not in set(order) and v not in (9770, 3183)]
    from confluent_eval5 import make_forward
    solve = make_forward(kind, info, seq, bestval)
    run, seq2 = make_slack_solver(kind, info, seq, bestval)

    # baseline forward-eval violation set
    base = solve(list(bestval), [])
    vbase = viol_atoms(A, base)
    print(f"forward-eval all-0 violated: {sorted(vbase)} ({len(vbase)})", flush=True)

    # try activating x_12779 via a single 22-side bit, then slack-active
    for bit in [None, 19520, 2795, 1858]:
        setb = [] if bit is None else [bit]
        v1 = solve(list(bestval), setb)  # pass 1
        x12779 = v1[12779]
        frozen = {24026: v1[18274]-v1[35186], 27116: v1[17728]-v1[1642]}
        v2 = run(list(v1), frozen)
        # after run, x_24026/x_27116 frozen; recheck the twist directly
        tw = (v2[9770]==v2[18274], v2[3183]==v2[17728])
        bad = viol_atoms(A, v2)
        print(f"\nbit={bit}: x_12779(pass1)={x12779}  frozen x_24026={frozen[24026]!=0} x_27116={frozen[27116]!=0}", flush=True)
        print(f"  slack-active twist: x_9770==x_18274 -> {tw[0]}, x_3183==x_17728 -> {tw[1]}", flush=True)
        print(f"  slack-active violated atoms: {len(bad)}: {sorted(bad)[:20]}", flush=True)
    print(f"\ndone ({time.time()-t0:.0f}s)", flush=True)

if __name__ == '__main__':
    main()
