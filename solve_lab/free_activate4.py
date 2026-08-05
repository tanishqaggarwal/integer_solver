#!/usr/bin/env python3
"""CLEAN mod-based free-var activation. On a setting with 642894|(x_18274-x_35186)
and 2|(x_17728-x_1642) [reachable ~16%], set:
  x_15=1; x_30077=(x_35186-x_18274)/642894; x_31807=(x_1642-x_17728)/2
  freeze x_24026=(x_18274-x_35186)/2; x_27116=(x_17728-x_1642)/2
Then x_9770=x_18274, x_3183=x_17728 EXACTLY, x_26977=0 (a1817 0=0), and a1813/a1815
SATISFIED (x_38215=x_30077, x_29437=x_31807). Only the ripple/verifier-square atoms
should break. Measure violations; try many mod-satisfying settings to find the
lowest. Verify any 0-violation state against ORIGINAL atoms -> witness."""
import json, time
from confluent_eval5 import build5, make_forward
from slack_active import make_slack_solver, viol_atoms
from propagate import NVARS
M=642894

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
    st=12345
    def rnd():
        nonlocal st; st=(st*6364136223846793005+1442695040888963407)&((1<<64)-1); return st>>33

    def activate(S, v0):
        x15=1
        x30077=(v0[35186]-v0[18274])//M
        x31807=(v0[1642]-v0[17728])//2
        frozen={24026:(v0[18274]-v0[35186])//2, 27116:(v0[17728]-v0[1642])//2}
        val=list(bestval); val[15]=x15; val[30077]=x30077; val[31807]=x31807
        v2=run(val, frozen)
        return v2

    best=None; ntried=0
    for _ in range(6000):
        if time.time()-t0>1500: break
        k=8+rnd()%30
        S=sorted(set(control[rnd()%len(control)] for _ in range(k)))
        v0=solve(list(bestval),S)
        if v0[12779]!=2 or v0[18274]==0: continue
        if (v0[18274]-v0[35186])%M!=0 or (v0[17728]-v0[1642])%2!=0: continue
        ntried+=1
        v2=activate(S, v0)
        tw=(v2[3183]==v2[17728], v2[9770]==v2[18274])
        bad=viol_atoms(A,v2)
        if best is None or len(bad)<best[0]:
            best=(len(bad), S, tw, sorted(bad))
            print(f"  mod-state #{ntried}: twist={tw}, x24026!=0={v2[24026]!=0}, violated {len(bad)}: {sorted(bad)[:16]} ({time.time()-t0:.0f}s)", flush=True)
            if len(bad)==0:
                allbad=viol_atoms(A,v2)  # original atoms (A is original)
                if not allbad:
                    json.dump({f"x_{i}":v2[i] for i in range(NVARS)}, open('cand_ACTIVATED_SOLVED.json','w'))
                    print("  *** SOLVED *** cand_ACTIVATED_SOLVED.json", flush=True); return
            if len(bad)<=10:
                json.dump({"bad":sorted(bad),"bits":S,"val":{str(i):v2[i] for i in range(NVARS)}}, open('activated_state.json','w'))
    print(f"\ntried {ntried} mod-satisfying states; best {best[0] if best else '?'} violated ({time.time()-t0:.0f}s)", flush=True)
    if best: print(f"  best twist={best[2]}, atoms={best[3][:16]}", flush=True)

if __name__=='__main__':
    main()
