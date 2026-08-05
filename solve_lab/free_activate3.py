#!/usr/bin/env python3
"""CLEAN free-var slack activation. Master free var x_15 grounds ALL slack cascades:
  x_37917=x_15 -> x_38215=x_30077*x_15 -> x_24026 -> x_3368 (9770 slack)
  x_7815 =x_15 -> x_29437=x_15*x_31807 -> x_27116 -> x_10466 (3183 slack)
  x_20510=x_15 -> x_26977=x_15*x_31302 (a1817 product slack for 9770 gap)
Set (x_12779=2 so x_14402=-1):
  x_15=1; x_31807=(x_1642-x_17728)/2  => x_3183=x_17728 (RIGID, needs mod-2)
  x_30077=knob                        => x_9770 = x_35186-642894*x_30077 (gap free)
  x_31302=6033033*(x_18274-x_9770)    => a1817 holds; Q40782 twist+a1817 terms CANCEL
So all a181x SATISFIED; remaining breakage = ripple / verifier squares R_i. Freeze
x_24026/x_27116 for ordering. Measure violations; scan x_30077 to shrink them."""
import json, time
from confluent_eval5 import build5, make_forward
from slack_active import make_slack_solver, viol_atoms
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
    run=make_slack_solver(kind,info,seq,bestval)[0]  # freezes {24026,27116}
    control=json.load(open('control_bits.json'))
    st=5
    def rnd():
        nonlocal st; st=(st*6364136223846793005+1442695040888963407)&((1<<64)-1); return st>>33

    # find x_12779=2 setting with (x_1642-x_17728) even and x_18274!=0
    S=None
    for _ in range(1500):
        k=8+rnd()%30
        cand=sorted(set(control[rnd()%len(control)] for _ in range(k)))
        v=solve(list(bestval),cand)
        if v[12779]==2 and v[18274]!=0 and (v[1642]-v[17728])%2==0:
            S=cand; v0=v; break
    if S is None:
        print("no suitable setting", flush=True); return
    print(f"setting k={len(S)}: x_12779=2, (x_1642-x_17728) even ({time.time()-t0:.0f}s)", flush=True)
    x31807=(v0[1642]-v0[17728])//2

    def activate(x30077):
        free={15:1, 30077:x30077, 31807:x31807, 31302:0}
        val=list(bestval)
        for f,x in free.items(): val[f]=x
        v1=solve(val, S)                              # pass1: get x_9770,x_18274,x_24026,x_27116
        x31302=6033033*(v1[18274]-v1[9770])           # a1817: x_26977=x_15*x_31302 (x_15=1)
        free[31302]=x31302
        val2=list(bestval)
        for f,x in free.items(): val2[f]=x
        v1b=solve(val2, S)
        frozen={24026:v1b[24026], 27116:v1b[27116], 26977:x31302}  # x_15=1 so x_26977=x_31302
        val3=list(bestval)
        for f,x in free.items(): val3[f]=x
        v2=run(val3, frozen)
        return v2

    best=None
    for x30077 in [1,2,-1,0,3,-2,7,-7,100,-100]:
        v2=activate(x30077)
        tw=(v2[3183]==v2[17728], v2[9770]==v2[18274])
        a1817ok = (v2[26977]==6033033*(v2[18274]-v2[9770]))
        bad=viol_atoms(A,v2)
        print(f"  x_30077={x30077}: x3183==x17728 {tw[0]}, a1817ok {a1817ok}, x24026!=0 {v2[24026]!=0}, violated {len(bad)}: {sorted(bad)[:14]}", flush=True)
        if best is None or len(bad)<best[0]: best=(len(bad), x30077, v2)
    print(f"\nbest: {best[0]} violated at x_30077={best[1]} ({time.time()-t0:.0f}s)", flush=True)
    if best[0]<=12:
        json.dump({f"x_{i}":best[2][i] for i in range(NVARS)}, open('free3_best.json','w'))
        print("saved free3_best.json", flush=True)

if __name__=='__main__':
    main()
