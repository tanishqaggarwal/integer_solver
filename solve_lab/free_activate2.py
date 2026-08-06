#!/usr/bin/env python3
"""CORE breakthrough attempt: activate BOTH slacks via the free roots x_15, x_30077,
x_31807 (both cascades ground at x_15). At x_12779=2 (x_14402=-1):
  x_38215 = x_30077*x_15,  x_24026 = -321447*x_38215,  x_3368 = 2*x_24026
  x_29437 = x_15*x_31807,  x_27116 = -x_29437,          x_10466 = 2*x_27116
Twist:  x_9770 = x_35186 - 642894*x_30077*x_15 = x_18274
        x_3183 = x_1642   - 2*x_15*x_31807     = x_17728
=> x_30077*x_15 = (x_35186-x_18274)/642894 =: P1 ;  x_15*x_31807 = (x_1642-x_17728)/2 =: P2.
Set x_15 = gcd(P1,P2), x_30077=P1/x_15, x_31807=P2/x_15. Cascade atoms stay SATISFIED
(a1813/a1815 hold since x_38215/x_29437 match). Freeze x_24026/x_27116 for ordering.
Report violations -- hope << 18."""
import json, time
from math import gcd
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

    def eval_free(bits, free, frozen=None):
        val=list(bestval)
        for f,x in free.items(): val[f]=x
        v1=solve(val, bits)
        if frozen is None: return v1
        val2=list(v1)
        for f,x in free.items(): val2[f]=x
        return run(val2, frozen)

    # search x_12779=2 settings; check P1,P2 integrality
    tried=0; good=0
    for _ in range(2500):
        k=8+rnd()%30
        S=sorted(set(control[rnd()%len(control)] for _ in range(k)))
        v=solve(list(bestval),S)
        if v[12779]!=2 or v[18274]==0: continue
        tried+=1
        d1=v[35186]-v[18274]; d2=v[1642]-v[17728]
        if d1%642894!=0 or d2%2!=0: continue
        P1=d1//642894; P2=d2//2
        good+=1
        g=gcd(P1,P2) or 1
        x15=g; x30077=P1//g if g else P1; x31807=P2//g if g else P2
        if x15==0: continue
        # pass1 to get activated x_24026/x_27116, then freeze
        v1=eval_free(S, {15:x15, 30077:x30077, 31807:x31807})
        frozen={24026:v1[24026], 27116:v1[27116]}
        v2=eval_free(S, {15:x15, 30077:x30077, 31807:x31807}, frozen)
        tw=(v2[9770]==v2[18274], v2[3183]==v2[17728])
        if tw[0] and tw[1]:
            bad=viol_atoms(A,v2)
            print(f"*** twist HOLDS via free-var activation: bits(k={len(S)}), x_15={x15}", flush=True)
            print(f"    a1813 satisfied? {24026 not in bad and 1813 not in bad}; violated atoms: {len(bad)}: {sorted(bad)[:20]}", flush=True)
            if len(bad)<=12:
                json.dump({"bad":sorted(bad),"bits":S,"free":{15:x15,30077:x30077,31807:x31807},
                           "val":{str(i):v2[i] for i in range(NVARS)}}, open('free_activate_state.json','w'))
                print("    saved state", flush=True)
            if len(bad)==0:
                json.dump({f"x_{i}":v2[i] for i in range(NVARS)}, open('cand_FREEACT_SOLVED.json','w'))
                print("    *** SOLVED ***", flush=True); return
            break
    print(f"searched: {tried} x_12779=2 states, {good} with P1,P2 integer ({time.time()-t0:.0f}s)", flush=True)

if __name__=='__main__':
    main()
