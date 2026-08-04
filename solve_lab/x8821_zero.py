#!/usr/bin/env python3
"""CREATIVE escape: the x_8821=0 regime. When x_8821=x_17810*x_27292=0, the div
wires collapse (a4954: x_6773=0, a13204: x_17233=0) and x_18274,x_17728 are FREED
from quantization -- determined by their LINEAR gates:
  x_18274 = x_25804 - x_35846   (a11387)   [or x_31434 - x_6283 (a11398)]
  x_17728 = x_27912 - x_28035   (a11388)
x_25804,x_35846,x_27912,x_28035 are UPSTREAM (independent of x_18274/x_17728) so
forward-eval computes them correctly even when x_8821=0. Read off the TRUE
x_18274/x_17728 there and test whether they can hit the 22-side x_9770/x_3183 values
(escaping the coprime trap). Also check x_6773==0, x_17233==0 (the regime's
constraints)."""
import json, time
from math import gcd
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
    b233=[b for b in control if b not in BITS22]; b22=[b for b in control if b in BITS22]

    base=solve(list(bestval),[])
    g=base[9770]; h=base[3183]
    # x_8821 is linear in ~18 bits; find bits that move it
    x8821movers=[b for b in b233 if solve(list(bestval),[b])[8821]!=base[8821]]
    print(f"x_8821 base={base[8821]}, moved by {len(x8821movers)} bits", flush=True)

    # true x_18274/x_17728 via linear gates (correct even at x_8821=0)
    def true18274(v): return v[25804]-v[35846]   # a11387
    def true17728(v): return v[27912]-v[28035]   # a11388

    st=71
    def rnd():
        nonlocal st; st=(st*6364136223846793005+1442695040888963407)&((1<<64)-1); return st>>33
    # search settings with x_8821=0 and measure true x_18274/x_17728 image
    found0=0; img18=set(); img17=set(); hit=[]
    # collect 22-side x_9770/x_3183 achievable values
    v9770=set(); v3183=set()
    for _ in range(120):
        k=1+rnd()%len(b22); S=sorted(set(b22[rnd()%len(b22)] for _ in range(k)))
        v=solve(list(bestval),S); v9770.add(v[9770]); v3183.add(v[3183])
    print(f"22-side: |x_9770 image|={len(v9770)}, |x_3183 image|={len(v3183)}", flush=True)

    for _ in range(400):
        k=1+rnd()%22; S=sorted(set(b233[rnd()%len(b233)] for _ in range(k)))
        v=solve(list(bestval),S)
        if v[8821]==0:
            found0+=1
            t18=true18274(v); t17=true17728(v)
            img18.add(t18); img17.add(t17)
            # regime constraints: x_6773 should be 0, x_17233 should be 0 (a4954/a13204 with x_8821=0)
            if t18 in v9770 and t18!=0: hit.append(('x18274=x9770', t18, S))
            if t17 in v3183 and t17!=0: hit.append(('x17728=x3183', t17, S))
    print(f"x_8821=0 settings found: {found0}/3000", flush=True)
    print(f"  true x_18274 image at x_8821=0: {len(img18)} distinct; sample {sorted(list(img18))[:4]}", flush=True)
    print(f"  true x_17728 image at x_8821=0: {len(img17)} distinct", flush=True)
    inter18=img18 & v9770; inter17=img17 & v3183
    print(f"  true x_18274(x_8821=0) ∩ x_9770(22) = {len(inter18)}: {sorted(inter18)[:4]}", flush=True)
    print(f"  true x_17728(x_8821=0) ∩ x_3183(22) = {len(inter17)}: {sorted(inter17)[:4]}", flush=True)
    # quantization check: do true values escape g2/h2?
    esc18=sum(1 for x in img18 if x!=0 and base[18274] and x%base[18274]!=0)
    esc17=sum(1 for x in img17 if x!=0 and base[17728] and x%base[17728]!=0)
    print(f"  true x_18274 escaping g2*Z: {esc18}/{len(img18)}; true x_17728 escaping h2*Z: {esc17}/{len(img17)}", flush=True)
    if hit: print(f"  *** TWIST HITS: {hit[:3]}", flush=True)
    print(f"done ({time.time()-t0:.0f}s)", flush=True)

if __name__=='__main__':
    main()
