#!/usr/bin/env python3
"""Q-5: structure of the 253 decoded leaf points inside the group."""
import sys, os, itertools, collections, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qgrp import *
N = 115792089237316195423570985008687907852837564279074904382605163141518161494337
L = leaves()
ks = sorted(L)
P = [L[k] for k in ks]
S = {pt:i for i,pt in enumerate(P)}
print('leaves', len(P), 'distinct', len(S))
print('target on curve:', oncur(TARGET), ' target is a leaf:', TARGET in S)
# embedding degree
for k in range(1,25):
    if pow(p,k,N)==1: print('embedding degree k =',k); break
else: print('embedding degree > 24')
# doubling ladder?
dbl = {}
for i,pt in enumerate(P):
    d = add(pt,pt)
    if d in S: dbl[i]=S[d]
print('leaves whose DOUBLE is also a leaf: %d/%d' % (len(dbl), len(P)))
# chains
succ = dbl
pred = collections.defaultdict(list)
for a,b in succ.items(): pred[b].append(a)
starts=[i for i in range(len(P)) if i not in pred]
chains=[]
for s in starts:
    c=[s]
    while c[-1] in succ: c.append(succ[c[-1]])
    chains.append(c)
chains.sort(key=len,reverse=True)
print('doubling chains (len):', [len(c) for c in chains][:12], 'total chains', len(chains))
# small-multiple relations between leaves
neg_ = {neg(pt):i for i,pt in enumerate(P)}
print('leaves whose negative is a leaf:', sum(1 for pt in P if pt in neg_ and neg_[pt]!=S[pt]))
json.dump({'ks':ks,'chains':chains}, open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'struct.json'),'w'))
