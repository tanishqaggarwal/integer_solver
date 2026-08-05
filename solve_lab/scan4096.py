#!/usr/bin/env python3
"""Scan all 2^12 combos of the 12 core-safe verifier selectors (x_2081=1 fixed).
For each: set selectors, load active pins to CONSTs, forward, close G1/G2, count fails.
The SAT agent only covered 6 of these -> ~4000 untested. Find any combo below agentA's 16."""
import heal_harness as H, json, itertools, sys
p=H.p
pins=json.load(open('pinrec.json'))
sels=[4287,5910,11368,12054,13195,16586,17406,18022,22562,23751,24365,28005]
pinmap={}
for pr in pins:
    ai,sel,tgt,cst,coef,hnd=pr[:6]
    if sel in sels:
        pinmap.setdefault(sel,{}).setdefault(tgt,cst)
pinmap={s:list(d.items()) for s,d in pinmap.items()}
vA=H.loadd('best_agentA_39022.json')
base={v:vA.get(v,0) for v in H.freeinp}

best=99; bestcombo=None; hits=[]
cnt=0
for combo in itertools.product([0,1],repeat=len(sels)):
    cnt+=1
    for v in H.freeinp: H.val[v]=base[v]
    H.val[2081]=1
    for s,bit in zip(sels,combo):
        H.val[s]=bit
        if bit==1:
            for tgt,cst in pinmap.get(s,[]):
                H.val[tgt]=cst%p
    H.forward()
    # close G1/G2
    H.val[17325]=0; H.val[9413]=0; H.forward()
    H.val[7068]=H.val[2099]; H.val[4432]=H.val[19964]; H.forward()
    nf=len(H.fails())
    if nf<best:
        best=nf; bestcombo=combo
    if nf<=15:
        hits.append((nf,combo))
    if cnt%512==0:
        print(f'  ...{cnt}/4096 scanned, best so far {best}', flush=True)
print(f'DONE. scanned {cnt}. BEST fails={best} combo={dict(zip(sels,bestcombo))}')
print(f'combos with <=15 fails: {len(hits)}')
for nf,combo in sorted(hits)[:20]:
    on=[s for s,b in zip(sels,combo) if b]
    print(f'  fails={nf}: selectors ON = {on}')
