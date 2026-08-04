#!/usr/bin/env python3
"""NEW: the div-by-zero that forces x_24026=0 happens ONLY at x_12779=1 (x_14402=0).
At x_12779=2,3,4 (x_14402=-1,-2,-3) forward-eval computes x_24026=321447*x_38215/
x_14402 CORRECTLY -> the 9770-side slack activates WITHIN forward-eval (no frozen
hack). x_38215=x_37917*x_30077 is controllable. Test: find bits giving x_12779=2,
measure x_24026/x_3368/x_9770, and whether x_9770 leaves g*Z (non-degenerate)."""
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
    base=solve(list(bestval),[])
    g=base[9770]
    print(f"base: x_12779={base[12779]}, x_24026={base[24026]}, x_9770={base[9770]}, x_38215={base[38215]}, x_14402={base[14402]}", flush=True)

    # find single bits and their effect on x_12779 and x_38215
    b12779={}; b38215={}
    for b in control:
        v=solve(list(bestval),[b])
        b12779[b]=v[12779]; b38215[b]=v[38215]
    # bits giving x_12779=2,3,4
    for target in (2,3,4):
        bits=[b for b in control if b12779[b]==target]
        print(f"single bits giving x_12779={target}: {bits[:12]} (n={len(bits)})", flush=True)
    b38=[b for b in control if b38215[b]!=base[38215]]
    print(f"single bits moving x_38215: {b38[:12]} (n={len(b38)})", flush=True)

    # try to get x_12779=2 AND x_38215!=0 via pairs/triples
    st=17
    def rnd():
        nonlocal st; st=(st*6364136223846793005+1442695040888963407)&((1<<64)-1); return st>>33
    got=[]
    b12779_2=[b for b in control if b12779[b]==2]
    print(f"\nsearching combos for x_12779 in {{2,3,4}} with x_24026 != 0 (slack active)...", flush=True)
    tried=0
    for _ in range(3000):
        k=1+rnd()%6
        S=sorted(set(control[rnd()%len(control)] for _ in range(k)))
        v=solve(list(bestval),S)
        tried+=1
        if v[12779] in (2,3,4) and v[24026]!=0:
            got.append((S, v[12779], v[24026], v[9770], v[3368], v[38215], v[14402]))
            if len(got)<=6:
                nd = "NONdegen" if (g==0 or v[9770]%g!=0) else "g-mult"
                print(f"  bits={S}: x_12779={v[12779]}, x_24026={v[24026]!=0}, x_3368={v[3368]!=0}, x_9770 {nd}", flush=True)
        if len(got)>=40: break
    print(f"\nfound {len(got)} slack-ACTIVE forward-eval states (x_12779 in 2..4, x_24026!=0)", flush=True)
    if got:
        # image of x_9770 in these states - does it escape g*Z / reach new values?
        vals9770=set(x[3] for x in got)
        esc=sum(1 for x in vals9770 if g and x%g!=0)
        print(f"  distinct x_9770 in slack-active states: {len(vals9770)}; escaping g*Z: {esc}", flush=True)
        print(f"  sample x_9770 values: {sorted(vals9770)[:4]}", flush=True)
    print(f"done ({time.time()-t0:.0f}s)", flush=True)

if __name__=='__main__':
    main()
