#!/usr/bin/env python3
"""CORE TEST: the escape cascade grounds at FREE vars x_15, x_30077 (df=None). Set
them nonzero -> forward-eval propagates x_37917=x_15, x_38215=x_30077*x_15,
x_24026=-321447*x_38215/x_14402 (needs x_14402 != 0, i.e. x_12779 != 1). If x_24026
activates, the slack is forward-eval-representable. Then tune x_15 so x_9770=x_18274.
Also find the 3183-side free roots (analog of x_15/x_30077 via x_29437)."""
import json, time
from confluent_eval5 import build5, make_forward
from slack_active import viol_atoms
from propagate import atom_vars, NVARS

def main():
    t0=time.time()
    A,kind,info,seq0,bestval,ncyc=build5()
    order=json.load(open('eval_order.json'))['order']
    defset=set(v for v in kind if kind[v]!='const')
    seq=[v for v in order if v in defset and v not in (9770,3183)]
    seq+=[v for v in (9770,3183) if v in defset]
    seq+=[v for v in defset if v not in set(order) and v not in (9770,3183)]
    solve=make_forward(kind,info,seq,bestval)
    control=json.load(open('control_bits.json'))

    # find a bit-setting with x_12779=2 (x_14402=-1) and x_18274 != 0
    st=5
    def rnd():
        nonlocal st; st=(st*6364136223846793005+1442695040888963407)&((1<<64)-1); return st>>33
    setbits=None
    for _ in range(800):
        k=8+rnd()%30
        S=sorted(set(control[rnd()%len(control)] for _ in range(k)))
        v=solve(list(bestval),S)
        if v[12779]==2 and v[18274]!=0 and v[35186]!=0:
            setbits=S; break
    if setbits is None:
        print("no suitable x_12779=2 bit-setting", flush=True); return
    v0=solve(list(bestval),setbits)
    print(f"x_12779={v0[12779]}, x_14402={v0[14402]}, x_18274!=0={v0[18274]!=0}, x_35186!=0={v0[35186]!=0}", flush=True)

    # Set free x_15, x_30077 and re-eval. First just x_15=1, x_30077=1 to see activation.
    def eval_with_free(bits, free):
        val=list(bestval)
        for f,x in free.items(): val[f]=x
        return solve(val, bits)
    v1=eval_with_free(setbits, {15:1, 30077:1})
    print(f"\nwith x_15=1, x_30077=1: x_37917={v1[37917]}, x_38215={v1[38215]}, x_24026={v1[24026]}, x_3368={v1[3368]}", flush=True)
    print(f"  x_9770 changed from {v0[9770]} to {v1[9770]}: {v0[9770]!=v1[9770]}", flush=True)

    # tune x_15 so x_9770 = x_18274:  x_9770 = x_35186 - 642894*x_30077*x_15  (x_12779=2)
    # => x_30077*x_15 = (x_35186 - x_18274)/642894 ; set x_30077=1, x_15 = that
    d=v0[35186]-v0[18274]
    print(f"\n(x_35186-x_18274)/642894 integer? {d%642894==0}", flush=True)
    if d%642894==0:
        x15=d//642894
        v2=eval_with_free(setbits, {15:x15, 30077:1})
        print(f"  set x_15={x15}: x_24026={v2[24026]!=0}, x_9770==x_18274? {v2[9770]==v2[18274]}", flush=True)
        bad=viol_atoms(A,v2)
        print(f"  violated atoms (twist may need 3183-side too): {len(bad)}: {sorted(bad)[:16]}", flush=True)
    print(f"done ({time.time()-t0:.0f}s)", flush=True)

if __name__=='__main__':
    main()
