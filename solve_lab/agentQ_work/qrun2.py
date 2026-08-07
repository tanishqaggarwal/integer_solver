#!/usr/bin/env python3
"""Q-10c: NON-CIRCULAR routing test.  Set ONLY the 256 selector bits.  The leaf coordinate wires
are *solved* from the pin atoms  sel*(w - C)  by the same unit propagation as everything else --
nothing is assigned by hand, so "OFF leaf behaves as the identity" is never presupposed.
Low weights (1,2,3,5,7) are tested explicitly, because that is the regime where the pass-through
branch is exercised at depth and where lowwt.py / wt7.py live."""
import json,random,sys,collections
sys.path.insert(0,'.')
import qsolve as Q
from qgrp import add,cs
p=Q.p
leaf={int(g):v for g,v in json.load(open('qleaf.json')).items()}
LP={g:(int(v[0]),int(v[1])) for g,v in leaf.items()}; LW={g:(v[2],v[3]) for g,v in leaf.items()}
lad=json.load(open('qladder.json')); e2s={int(k):v for k,v in lad['exp2sel'].items()}
ST=[x for x in json.load(open('qstages.json'))['stages'] if 'u3' in x]
ROOTX,ROOTY=24468,18956
G=[json.loads(l) for l in open('atoms/gates.jsonl')]
tg={d['t'] for d in G}; free=[v for v in range(Q.NV) if v not in tg]
def fold(S):
    F=None
    for g in S: F=add(F,LP[g])
    return F
print('%-5s | %-11s %-11s | %-11s %-11s | %-9s %-7s'%(
      'w','ON X solved','ON X correct','OFF X solved','OFF Y = 0','gadgets','root'))
random.seed(7); sel=sorted(LP)
for w in (1,2,3,5,7,128):
    S=set(random.sample(sel,w))
    V=[None]*Q.NV
    for g in LP: V[g]=1 if g in S else 0          # ONLY the selectors are set
    s,c=Q.propagate(V)
    onx=[g for g in S if V[LW[g][0]] is not None]
    onok=[g for g in onx if V[LW[g][0]]==(LP[g][0]-cs)%p and V[LW[g][1]]==LP[g][1]%p]
    off=[g for g in LP if g not in S]
    offx=[g for g in off if V[LW[g][0]] is not None]
    offy0=[g for g in off if V[LW[g][1]]==0]
    ng=sum(1 for x in ST if V[x['u3']] is not None)
    print('%-5d | %-11s %-11s | %-11s %-11s | %-9s %-7s'%(
        w,'%d/%d'%(len(onx),len(S)),'%d/%d'%(len(onok),len(S)),
        '%d/%d'%(len(offx),len(off)),'%d/%d'%(len(offy0),len(off)),
        '%d/383'%ng, V[ROOTX] is not None))
    if w==1:
        print('      free inputs solved: %d/%d ; wires known %d/%d ; contradictions %d'%(
            sum(1 for v in free if V[v] is not None),len(free),
            sum(1 for v in V if v is not None),Q.NV,c))
