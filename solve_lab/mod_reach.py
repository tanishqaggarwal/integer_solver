#!/usr/bin/env python3
"""Is the modular twist constraint reachable? Need x_18274 ≡ x_35186 (mod 642894)
[9770 side] and x_17728 ≡ x_1642 (mod 2) [3183 side], with x_12779=2. Measure the
image of (x_18274-x_35186) mod 642894 and mod its prime factors, and whether 0 is
reachable. Small modulus (~20 bits) => if controllable, SOLVABLE by targeted search /
CRT, unlike the full 296-bit collision."""
import json, time
from math import gcd
from confluent_eval5 import build5, make_forward
M=642894  # = 2*3*7*15307

def factor(n):
    f={}; d=2
    while d*d<=n:
        while n%d==0: f[d]=f.get(d,0)+1; n//=d
        d+=1
    if n>1: f[n]=f.get(n,0)+1
    return f

def main():
    t0=time.time()
    print(f"642894 = {factor(M)}", flush=True)
    A,kind,info,seq0,bestval,ncyc=build5()
    order=json.load(open('eval_order.json'))['order']
    defset=set(v for v in kind if kind[v]!='const')
    seq=[v for v in order if v in defset and v not in (9770,3183)]
    seq+=[v for v in (9770,3183) if v in defset]
    seq+=[v for v in defset if v not in set(order) and v not in (9770,3183)]
    solve=make_forward(kind,info,seq,bestval)
    control=json.load(open('control_bits.json'))
    st=5
    def rnd():
        nonlocal st; st=(st*6364136223846793005+1442695040888963407)&((1<<64)-1); return st>>33
    base=solve(list(bestval),[])
    print(f"base (x_18274-x_35186) mod 642894 = {(base[18274]-base[35186])%M}", flush=True)

    img=set(); imgmod2=set(); reach0=0; reach0both=0; n2=0
    samples=[]
    for _ in range(2500):
        k=8+rnd()%30
        S=sorted(set(control[rnd()%len(control)] for _ in range(k)))
        v=solve(list(bestval),S)
        if v[12779]!=2: continue
        n2+=1
        r=(v[18274]-v[35186])%M
        r2=(v[17728]-v[1642])%2
        img.add(r); imgmod2.add(r2)
        samples.append((S,r,r2))
        if r==0: reach0+=1
        if r==0 and r2==0: reach0both+=1
    print(f"x_12779=2 states: {n2}", flush=True)
    print(f"|image of (x_18274-x_35186) mod 642894| = {len(img)}; contains 0? {0 in img}", flush=True)
    print(f"  sample residues: {sorted(img)[:12]}", flush=True)
    print(f"(x_17728-x_1642) mod 2 image: {sorted(imgmod2)}", flush=True)
    print(f"states with mod-642894 ≡ 0: {reach0}; with BOTH constraints: {reach0both}", flush=True)
    # per-prime-factor reachability
    for p in factor(M):
        imgp=set((v% p) for v in img) if False else set()
        # recompute residues mod p from samples
        rp=set((r % p) for _,r,_ in samples)
        print(f"  mod {p}: image {sorted(rp)[:10]} (n={len(rp)}); contains 0? {0 in rp}", flush=True)
    print(f"done ({time.time()-t0:.0f}s)", flush=True)

if __name__=='__main__':
    main()
