#!/usr/bin/env python3
"""CREATIVE: what is the polynomial DEGREE of the numerator x_6773 (=g2*x18274*x8821
region) and x_17233 in the 233 bits? Pure-linear fits only 5-7/50. Test the
DEGREE-2 model x(S) = base + sum_i d_i + sum_{i<j in S} d_ij by measuring pairwise
interactions d_ij on demand for each test subset. If degree-2 fits, the twist
collision x_18274 = x_9770 becomes a QUADRATIC Diophantine (tractable), and the d_ij
that ESCAPE the quantum g2 are the mechanism that breaks coprime quantization."""
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
    b233=[b for b in control if b not in BITS22]
    base=solve(list(bestval),[])
    g2=base[18274]

    for W in (6773, 17233, 18274, 17728):
        wbase=base[W]
        d={}  # single deltas
        for b in b233:
            d[b]=solve(list(bestval),[b])[W]-wbase
        movers=[b for b in b233 if d[b]!=0]
        # pairwise interaction cache
        dij={}
        def inter(i,j):
            key=(min(i,j),max(i,j))
            if key not in dij:
                v=solve(list(bestval),[i,j])[W]
                dij[key]=v-wbase-d[i]-d[j]
            return dij[key]
        # test degree-1 and degree-2 on random subsets of movers
        st=20+W
        def rnd():
            nonlocal st; st=(st*6364136223846793005+1442695040888963407)&((1<<64)-1); return st>>33
        lin_ok=deg2_ok=ntest=0
        esc_inter=0; inter_seen=0
        for _ in range(30):
            k=2+rnd()%6
            S=sorted(set(movers[rnd()%len(movers)] for _ in range(k)))
            if len(S)<2: continue
            ntest+=1
            actual=solve(list(bestval),S)[W]
            lin_pred=wbase+sum(d[b] for b in S)
            if actual==lin_pred: lin_ok+=1
            # degree-2 pred
            q=wbase+sum(d[b] for b in S)
            for a in range(len(S)):
                for b in range(a+1,len(S)):
                    ij=inter(S[a],S[b]); q+=ij
                    inter_seen+=1
                    if ij!=0 and (g2==0 or ij%g2!=0): esc_inter+=1
            if actual==q: deg2_ok+=1
        print(f"x_{W}: {len(movers)} movers; linear fit {lin_ok}/{ntest}; degree-2 fit {deg2_ok}/{ntest}", flush=True)
        nzint=[v for v in dij.values() if v!=0]
        if nzint:
            gi=0
            for v in nzint: gi=gcd(gi,abs(v))
            print(f"   nonzero pair-interactions: {len(nzint)}/{len(dij)}; gcd={gi} (bits {gi.bit_length()}); gcd(inter,g2)={gcd(gi,g2)}", flush=True)
            print(f"   interactions ESCAPING g2*Z: {esc_inter}/{inter_seen}", flush=True)
        print(f"   ({time.time()-t0:.0f}s)", flush=True)
    print(f"done ({time.time()-t0:.0f}s)", flush=True)

if __name__=='__main__':
    main()
