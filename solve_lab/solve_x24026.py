#!/usr/bin/env python3
"""NEW: x_24026 is a FREE continuous knob (the 9770 side is bridged by a1817's
product slack x_26977, not a rigid constraint). Instead of forcing x_9770=x_18274,
SOLVE the verifier squares for x_24026. a1817 is auto-satisfied (x_26977 :=
6033033*(x_18274-x_9770)). Only a44271 (x_3183=x_17728) is rigid -> x_27116 fixed.

For a slack-active state parameterized by x_24026 (x_27116 held at x_17728-x_1642),
test whether Q40782 and the other broken atoms are LINEAR (or low-degree) in
x_24026; if so solve Q40782(x_24026)=0 exactly and re-check all atoms."""
import json, time
from confluent_eval5 import build5, make_forward
from slack_active import make_slack_solver, viol_atoms
from check_square import try_sqrt
from propagate import NVARS

def main():
    t0=time.time()
    A,kind,info,seq0,bestval,ncyc=build5()
    order=json.load(open('eval_order.json'))['order']
    defset=set(v for v in kind if kind[v]!='const')
    seq=[v for v in order if v in defset and v not in (9770,3183)]
    seq+=[v for v in (9770,3183) if v in defset]
    seq+=[v for v in defset if v not in set(order) and v not in (9770,3183)]
    solve=make_forward(kind,info,seq,bestval)
    run,seq2=make_slack_solver(kind,info,seq,bestval)
    Q40782=try_sqrt(A[40782])
    def qval(val):
        s=0
        for m,c in Q40782.items():
            t=c
            for x in m: t*=val[x]
            s+=t
        return s

    ACT=1858
    v1=solve(list(bestval),[ACT])
    x27116=v1[17728]-v1[1642]      # rigid: x_3183=x_17728
    x9770_target=v1[18274]         # baseline (x_24026 = x18274-x35186 => x9770=x18274)
    base24026=v1[18274]-v1[35186]

    # evaluate Q40782 as a function of x_24026 at several points
    def state(x24026):
        return run(list(v1), {24026:x24026, 27116:x27116})
    pts=[]
    for d in (0, 1, 2, 1000, base24026//2 if base24026 else 7):
        x=base24026+d if d<1000 else d
        val=state(x)
        pts.append((x, qval(val), val[9770]==val[18274], val[3183]==val[17728]))
    print("x_24026 offset -> Q40782 residual (twist9770,twist3183):", flush=True)
    for x,q,t9,t3 in pts:
        print(f"   x24026={'base' if x==base24026 else x}: Q={q}  (tw9={t9}, tw3={t3})", flush=True)
    # linearity in x_24026: Q(base)-Q(base+1) etc
    q0=pts[0][1]; q1=pts[1][1]; q2=pts[2][1]
    slope=q1-q0; curv=(q2-q1)-(q1-q0)
    print(f"\nQ slope (per +1 x_24026) = {slope}", flush=True)
    print(f"Q curvature (2nd diff) = {curv}  ({'LINEAR' if curv==0 else 'nonlinear'})", flush=True)
    if curv==0 and slope!=0:
        # solve Q(base + k) = 0 : q0 + k*slope = 0 -> k = -q0/slope
        if (-q0)%slope==0:
            k=(-q0)//slope
            xsol=base24026+k
            val=state(xsol)
            print(f"SOLVED Q40782=0 at x_24026=base+{k}; Q={qval(val)}", flush=True)
            bad=viol_atoms(A,val)
            print(f"  total violated atoms now: {len(bad)}: {sorted(bad)[:20]}", flush=True)
            if not bad:
                json.dump({f'x_{i}':val[i] for i in range(NVARS)},open('cand_x24026_SOLVED.json','w'))
                print("  *** SOLVED ***", flush=True)
        else:
            print(f"  Q40782=0 has no integer x_24026 (q0={q0} not divisible by slope={slope})", flush=True)
    print(f"done ({time.time()-t0:.0f}s)", flush=True)

if __name__=='__main__':
    main()
