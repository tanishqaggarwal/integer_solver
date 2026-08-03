#!/usr/bin/env python3
"""Is x_17233 (= x_8821*x_17728, the numerator) LINEAR in the 233 bits? If so and
x_8821 can be pinned to 1, then x_17728 = x_17233 = base + sum b_i*r_i, and
'x_17728 = target' is an integer subset-sum -> LLL-attackable if density<0.94.
Extract residues r_i, test linearity, compute density. Same for x_6773 (x_18274
numerator). Also: how many 233-bits control x_8821 vs the numerator (can we pin
x_8821=1 and still steer x_17233 over a large range)?"""
import json, time
from math import gcd
from functools import reduce
from confluent_eval5 import build5, make_forward
BITS22=set([1782,1858,2795,2800,3483,5443,10652,19520,21188,21588,23634,26947,27512,29682,30104,30596,30658,30792,33251,37748,37885,38116])

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
    b233=[b for b in control if b not in BITS22]
    base=solve(list(bestval),[])

    for W,name in [(17233,'x_17233 (num of x_17728)'),(6773,'x_6773 (num of x_18274)'),(8821,'x_8821 (denom)')]:
        delta={}
        for b in b233:
            delta[b]=solve(list(bestval),[b])[W]-base[W]
        st=55
        def rnd():
            nonlocal st; st=(st*6364136223846793005+1442695040888963407)&((1<<64)-1); return st>>33
        nlin=0;ntest=50
        for _ in range(ntest):
            k=1+rnd()%14
            S=sorted(set(b233[rnd()%len(b233)] for _ in range(k)))
            v=solve(list(bestval),S)
            pred=base[W]+sum(delta[b] for b in S)
            if v[W]==pred: nlin+=1
        nz=[d for d in delta.values() if d!=0]
        movers=[b for b in b233 if delta[b]!=0]
        g=reduce(gcd,[abs(d) for d in nz]) if nz else 0
        print(f"\n{name}: linear {nlin}/{ntest}; {len(nz)} bits move it; base bits={base[W].bit_length()}", flush=True)
        if nz:
            red=[abs(d)//g for d in nz]
            maxb=max(r.bit_length() for r in red);
            print(f"   gcd(deltas) bits={g.bit_length()}; deltas/gcd max bits={maxb}; density(n/maxbits)={len(nz)/max(maxb,1):.2f}", flush=True)
        print(f"   movers: {movers[:8]}... (n={len(movers)})", flush=True)
    print(f"\ndone ({time.time()-t0:.0f}s)", flush=True)

if __name__=='__main__':
    main()
