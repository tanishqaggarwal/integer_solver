#!/usr/bin/env python3
"""Target the DEGENERATE-ZERO twist: x_3183=x_17728=0 AND x_9770=x_18274=0, which
satisfies all twist atoms with NO slack ripple. Search 22-side bits zeroing
(x_3183,x_9770) and 233-side bits zeroing (x_17728,x_18274); combine and count
total violated atoms. If low, the degenerate witness is near-reachable (repairable);
this avoids the huge-slack ripple entirely."""
import json, time
from confluent_eval5 import build5, make_forward
BITS22=set([1782,1858,2795,2800,3483,5443,10652,19520,21188,21588,23634,26947,27512,29682,30104,30596,30658,30792,33251,37748,37885,38116])

def viol(A, val):
    bad=[]
    for a,poly in enumerate(A):
        s=0
        for m,c in poly.items():
            t=c
            for x in m: t*=val[x]
            s+=t
        if s: bad.append(a)
    return bad

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
    b22=[b for b in control if b in BITS22]
    b233=[b for b in control if b not in BITS22]
    st=101
    def rnd():
        nonlocal st; st=(st*6364136223846793005+1442695040888963407)&((1<<64)-1); return st>>33

    # search 22-side for x_3183=0 and/or x_9770=0
    best22=None  # (setbits, x3183, x9770)
    for _ in range(4000):
        k=1+rnd()%len(b22)
        S=sorted(set(b22[rnd()%len(b22)] for _ in range(k)))
        v=solve(list(bestval),S)
        score=(v[3183]==0)+(v[9770]==0)
        if best22 is None or score>best22[0] or (score==best22[0] and abs(v[3183])+abs(v[9770])<best22[3]):
            best22=(score,S,(v[3183],v[9770]),abs(v[3183])+abs(v[9770]))
    print(f"22-side best: score={best22[0]} bits={best22[1]} (x3183,x9770)={best22[2]}", flush=True)
    # search 233-side for x_17728=0 and/or x_18274=0
    best233=None
    for _ in range(6000):
        k=1+rnd()%20
        S=sorted(set(b233[rnd()%len(b233)] for _ in range(k)))
        v=solve(list(bestval),S)
        score=(v[17728]==0)+(v[18274]==0)
        if best233 is None or score>best233[0] or (score==best233[0] and abs(v[17728])+abs(v[18274])<best233[3]):
            best233=(score,S,(v[17728],v[18274]),abs(v[17728])+abs(v[18274]))
    print(f"233-side best: score={best233[0]} bits={best233[1][:10]}... (x17728,x18274)={best233[2]}", flush=True)

    # combine
    S=sorted(set(best22[1])|set(best233[1]))
    v=solve(list(bestval),S)
    print(f"\ncombined: x3183={v[3183]}, x9770={v[9770]}, x17728={v[17728]}, x18274={v[18274]}", flush=True)
    tw=(v[3183]==v[17728], v[9770]==v[18274])
    print(f"  twist holds: {tw}", flush=True)
    bad=viol(A,v)
    print(f"  total violated atoms: {len(bad)}: {sorted(bad)[:15]}", flush=True)
    print(f"done ({time.time()-t0:.0f}s)", flush=True)

if __name__=='__main__':
    main()
