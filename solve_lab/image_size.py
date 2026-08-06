#!/usr/bin/env python3
"""CREATIVE: is the 233-side image of x_17728 (and x_18274) SMALL? If it saturates
at a few hundred distinct values (like the 22-side's ~45 for x_3183), enumerate both
sides and find the twist collision DIRECTLY -- no trapdoor inversion. Sample with
HIGH bit-weights (not just low-weight) to probe the true image. Track saturation."""
import json, time
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
    b22=[b for b in control if b in BITS22]; b233=[b for b in control if b not in BITS22]
    st=424242
    def rnd():
        nonlocal st; st=(st*6364136223846793005+1442695040888963407)&((1<<64)-1); return st>>33

    # 22-side images (should saturate ~45/27)
    i3183=set(); i9770=set()
    for n in range(1200):
        w=rnd()%(len(b22)+1); S=[b for b in b22 if rnd()&1] if w==0 else sorted(set(b22[rnd()%len(b22)] for _ in range(w)))
        v=solve(list(bestval),S); i3183.add(v[3183]); i9770.add(v[9770])
    print(f"22-side after 1200 samples: |x_3183|={len(i3183)}, |x_9770|={len(i9770)}", flush=True)

    # 233-side images with HIGH weight sampling; track saturation
    i17728=set(); i18274=set()
    checkpts=[500,1000,2000,3500]
    ci=0
    for n in range(1,3501):
        w=1+rnd()%len(b233)               # weight 1..233 (often high)
        S=sorted(set(b233[rnd()%len(b233)] for _ in range(w)))
        v=solve(list(bestval),S); i17728.add(v[17728]); i18274.add(v[18274])
        if ci<len(checkpts) and n==checkpts[ci]:
            print(f"  233-side @ {n} samples: |x_17728|={len(i17728)}, |x_18274|={len(i18274)} ({time.time()-t0:.0f}s)", flush=True)
            ci+=1
    # collision check
    c1=i3183 & i17728; c2=i9770 & i18274
    print(f"\nCOLLISIONS: x_3183∩x_17728 = {len(c1)} {sorted(x for x in c1)[:5]}", flush=True)
    print(f"           x_9770∩x_18274 = {len(c2)} {sorted(x for x in c2)[:5]}", flush=True)
    nz1=[x for x in c1 if x!=0]; nz2=[x for x in c2 if x!=0]
    print(f"  NONZERO collisions: twist3183={len(nz1)}, twist9770={len(nz2)}", flush=True)
    if nz1 or nz2: print(f"  *** NONZERO TWIST COLLISION FOUND: {nz1[:3]} {nz2[:3]}", flush=True)
    print(f"done ({time.time()-t0:.0f}s)", flush=True)

if __name__=='__main__':
    main()
