#!/usr/bin/env python3
"""Alternative escape: make the twist hold from the 233-SIDE. Freeze x_18274 :=
x_9770(native 22-side) and x_17728 := x_3183(native), instead of activating the
22-side slacks. The ripple then lands on the 233-side (x_6773/x_17233/x_8821 web)
rather than the 22-side slack web -- possibly breaking different/fewer verifier
squares. Compare violation count to the slack-active approach (18)."""
import json, time
from confluent_eval5 import build5, make_forward
from slack_active import viol_atoms
from propagate import atom_vars, NVARS

def make_frz_solver(kind, info, seq, bestval, FREEZE):
    seq2 = [v for v in seq if v not in FREEZE]
    def run(val, frozen):
        for v in FREEZE: val[v] = frozen.get(v, val[v])
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
    return run

def main():
    t0=time.time()
    A,kind,info,seq0,bestval,ncyc=build5()
    order=json.load(open('eval_order.json'))['order']
    defset=set(v for v in kind if kind[v]!='const')
    seq=[v for v in order if v in defset and v not in (9770,3183)]
    seq+=[v for v in (9770,3183) if v in defset]
    seq+=[v for v in defset if v not in set(order) and v not in (9770,3183)]
    solve=make_forward(kind,info,seq,bestval)
    # freeze x_18274, x_17728 (233-side twist activation)
    run=make_frz_solver(kind,info,seq,bestval,{18274,17728})
    control=json.load(open('control_bits.json'))

    for bits in [[], [1858]]:
        v1=solve(list(bestval),bits)
        frozen={18274:v1[9770], 17728:v1[3183]}
        v2=run(list(v1),frozen)
        tw=(v2[9770]==v2[18274], v2[3183]==v2[17728])
        bad=viol_atoms(A,v2)
        print(f"bits={bits}: freeze x_18274:={v1[9770]!=0}, x_17728:={v1[3183]!=0}", flush=True)
        print(f"  233-side twist: x9770==x18274 {tw[0]}, x3183==x17728 {tw[1]}", flush=True)
        print(f"  violated atoms: {len(bad)}: {sorted(bad)[:16]}", flush=True)
    print(f"done ({time.time()-t0:.0f}s)", flush=True)

if __name__=='__main__':
    main()
