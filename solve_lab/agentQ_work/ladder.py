#!/usr/bin/env python3
"""Q-6: link the doubling chains into a single ladder and locate the target's index if easy."""
import sys, os, json, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qgrp import *
N = 115792089237316195423570985008687907852837564279074904382605163141518161494337
Lv = leaves(); ks=sorted(Lv); P=[Lv[k] for k in ks]; S={pt:i for i,pt in enumerate(P)}
succ={}
for i,pt in enumerate(P):
    d=add(pt,pt)
    if d in S: succ[i]=S[d]
pred=collections.defaultdict(list)
for a,b in succ.items(): pred[b].append(a)
starts=[i for i in range(len(P)) if not pred[i]]
chains=[]
for s in starts:
    c=[s]
    while c[-1] in succ: c.append(succ[c[-1]])
    chains.append(c)
chains.sort(key=len,reverse=True)
print('chains:', [len(c) for c in chains])
heads={c[0]:ci for ci,c in enumerate(chains)}
# try to join tail->head with 1,2,3 missing doublings
for ci,c in enumerate(chains):
    t=P[c[-1]]
    for gap in range(1,5):
        q=mul(pow(2,gap),t)
        if q in S and S[q] in heads:
            print('chain %d (len %d) --2^%d--> head of chain %d (len %d)'%(ci,len(c),gap,heads[S[q]],len(chains[heads[S[q]]])))
# build full ladder: order chains
order=[]
rem=list(range(len(chains)))
link={}
for ci,c in enumerate(chains):
    t=P[c[-1]]
    for gap in range(1,5):
        q=mul(pow(2,gap),t)
        if q in S and S[q] in heads: link[ci]=(heads[S[q]],gap); break
print('links:', link)
tails=set(link.keys()); tgts={v[0] for v in link.values()}
first=[ci for ci in range(len(chains)) if ci not in tgts]
print('chain with no predecessor:', first)
if len(first)==1 and len(link)==len(chains)-1:
    seq=[first[0]]
    while seq[-1] in link: seq.append(link[seq[-1]][0])
    print('chain order:', seq, 'lens', [len(chains[i]) for i in seq])
    # exponents
    G=P[chains[seq[0]][0]]
    e=0; ladder={}
    for si,ci in enumerate(seq):
        for j,idx in enumerate(chains[ci]):
            ladder[e]=idx; e+=1
        if ci in link: e += link[ci][1]-1
    print('ladder spans exponents 0..%d, %d present, %d missing'%(e-1,len(ladder),e-len(ladder)))
    miss=[i for i in range(e) if i not in ladder]
    print('missing exponents:', miss)
    # verify
    ok=all(mul(pow(2,i),G)==P[ladder[i]] for i in sorted(ladder))
    print('ladder verified L_i == 2^i * G :', ok)
    json.dump({'G_leafvar':ks[chains[seq[0]][0]],
               'ladder':{str(i):ks[ladder[i]] for i in ladder},
               'missing':miss,'N':str(N)},
              open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'ladder.json'),'w'))
    print('G leaf var =', ks[chains[seq[0]][0]], 'G =', G)
    print('TARGET =', TARGET)
