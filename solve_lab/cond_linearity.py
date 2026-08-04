#!/usr/bin/env python3
"""Is x_17233 (numerator of x_17728) a clean SUBSET-SUM over the 193 numerator bits
when the x_8821 (denominator) bits are held fixed? The unconditional nonlinearity
(7/50 linear) may come from x_8821 varying. If conditionally linear, x_17233 = base +
sum b_i*r_i is a subset-sum with density ~193/296<0.94 -> LLL-attackable. Identify
the x_8821 bits, hold them at 0, and test linearity over the rest."""
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
    # x_8821 bits = those that move x_8821
    x8821bits=[b for b in b233 if solve(list(bestval),[b])[8821]!=base[8821]]
    numbits=[b for b in b233 if b not in x8821bits]
    print(f"x_8821 bits: {len(x8821bits)}; numerator bits: {len(numbits)}", flush=True)

    # deltas of x_17233 over numbits (x_8821 bits held 0)
    W=17233
    delta={}
    for b in numbits:
        delta[b]=solve(list(bestval),[b])[W]-base[W]
    st=71
    def rnd():
        nonlocal st; st=(st*6364136223846793005+1442695040888963407)&((1<<64)-1); return st>>33
    lin=0; ntest=40
    for _ in range(ntest):
        k=2+rnd()%10
        S=sorted(set(numbits[rnd()%len(numbits)] for _ in range(k)))
        v=solve(list(bestval),S)  # only numerator bits set, x_8821 bits stay 0
        pred=base[W]+sum(delta[b] for b in S)
        if v[W]==pred: lin+=1
    nz=[d for d in delta.values() if d!=0]
    print(f"x_17233 CONDITIONAL linearity (x_8821 bits=0): {lin}/{ntest}", flush=True)
    if nz:
        g=reduce(gcd,[abs(d) for d in nz])
        red=[abs(d)//g for d in nz]
        print(f"  {len(nz)} moving bits; gcd(deltas) bits={g.bit_length()}; deltas/gcd max bits={max(r.bit_length() for r in red)}", flush=True)
        print(f"  if linear: subset-sum density = {len(nz)}/{base[W].bit_length()} = {len(nz)/base[W].bit_length():.2f}", flush=True)
    # also test x_6773 (num of x_18274) conditional
    W2=6773; d2={}
    for b in numbits: d2[b]=solve(list(bestval),[b])[W2]-base[W2]
    lin2=0
    for _ in range(ntest):
        k=2+rnd()%10
        S=sorted(set(numbits[rnd()%len(numbits)] for _ in range(k)))
        v=solve(list(bestval),S)
        if v[W2]==base[W2]+sum(d2[b] for b in S): lin2+=1
    print(f"x_6773 CONDITIONAL linearity (x_8821 bits=0): {lin2}/{ntest}", flush=True)
    print(f"done ({time.time()-t0:.0f}s)", flush=True)

if __name__=='__main__':
    main()
