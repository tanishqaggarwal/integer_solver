#!/usr/bin/env python3
"""Q-9b: with all 256 leaves decoded straight from EQUATIONS.txt, re-verify the doubling ladder."""
import json,sys,os
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from qgrp import add,neg,mul,oncur,p
L={int(g):(int(v[0]),int(v[1])) for g,v in json.load(open('qleaf.json')).items()}
S={pt:g for g,pt in L.items()}
print('leaves %d  distinct points %d  all on curve %s'%(len(L),len(S),all(oncur(pt) for pt in L.values())))
dbl={g:S[add(pt,pt)] for g,pt in L.items() if add(pt,pt) in S}
print('leaves whose double is also a leaf: %d/256'%len(dbl))
pred={}
for a,b in dbl.items(): pred.setdefault(b,[]).append(a)
starts=[g for g in L if g not in pred]
chains=[]
for s in starts:
    c=[s]
    while c[-1] in dbl: c.append(dbl[c[-1]])
    chains.append(c)
chains.sort(key=len,reverse=True)
print('doubling chains:',[len(c) for c in chains])
# stitch: each chain tail -> next chain head by one extra doubling
order=[]
rem=chains[:]
cur=rem.pop(0); order=cur[:]
while rem:
    t=L[order[-1]]; nxt=add(t,t)
    hit=None
    for ci,c in enumerate(rem):
        if L[c[0]]==nxt: hit=ci; break
    if hit is None:
        # allow one *missing* leaf between: order[-1] doubled twice
        n2=add(nxt,nxt)
        for ci,c in enumerate(rem):
            if L[c[0]]==n2: hit=ci; break
        if hit is None: print('STITCH FAILED at',order[-1]); break
        order.append(None)
    order += rem.pop(hit)
print('stitched length (None = leaf absent):',len(order),' Nones:',order.count(None))
G=L[order[0]]
ok=sum(1 for i,g in enumerate(order) if g is not None and L[g]==mul(pow(2,i),G))
print('leaves matching 2^i * G for i = 0..%d : %d / %d'%(len(order)-1,ok,len(order)-order.count(None)))
exps={}
for i,g in enumerate(order):
    if g is not None: exps[i]=g
json.dump({'G_leafvar':order[0],'exp2sel':{str(i):g for i,g in exps.items()},
           'sel2exp':{str(g):i for i,g in exps.items()}},open('qladder.json','w'))
print('G selector =',order[0])
