#!/usr/bin/env python3
"""Q-10d: the degenerate branch.  Every gadget was verified for DISTINCT input points.  Here the two
inputs are set EQUAL (P_a = P_b), where the chord is undefined and the group law needs the tangent.
If the gadget cannot represent doubling, then two equal live inputs are simply infeasible."""
import json,collections,random,sys
sys.path.insert(0,'.')
from qgrp import add,neg,p,cs,A_,B_
import re
G=[json.loads(l) for l in open('atoms/gates.jsonl')]
defs={}
for d in G: defs.setdefault(d['t'],d)
ST=[x for x in json.load(open('qstages.json'))['stages'] if 'u3' in x]
def cone(root,cut):
    seen=set(); stk=[root]
    while stk:
        w=stk.pop()
        if w in seen or w in cut: continue
        seen.add(w); d=defs.get(w)
        if d:
            for v in d['vids']: stk.append(v)
    done=set(cut); out=[]; prog=True
    while prog:
        prog=False
        for w in seen:
            if w in done or w not in defs: continue
            if all(v in done for v in defs[w]['vids']): out.append(w); done.add(w); prog=True
    return out
def ev(order,val):
    for w in order:
        d=defs[w]; e=d['rhs']
        for v in sorted(set(d['vids']),reverse=True): e=e.replace('x_%d'%v,'(%d)'%val[v])
        val[w]=eval(e,{'__builtins__':{}})%p
random.seed(2)
def randpt():
    while True:
        x=random.randrange(p); r=(pow(x,3,p)+A_*x+B_)%p
        y=pow(r,(p+1)//4,p)
        if y*y%p==r: return (x,y)
ok=0; nz=0; tested=0
for s in ST:
    cut={s['ua'],s['ub'],s['ya'],s['yb']}
    o1=cone(s['R1'],cut)
    P=randpt(); D=add(P,P)
    val=collections.defaultdict(int)
    val[s['ua']]=val[s['ub']]=(P[0]-cs)%p
    val[s['ya']]=val[s['yb']]=P[1]%p
    val[s['u3']]=(D[0]-cs)%p; val[s['y3']]=D[1]%p
    ev(o1,val); tested+=1
    if val[s['R1']]%p==0: ok+=1
    else: nz+=1
print('gadgets fed two EQUAL live points, output set to the true double 2P:')
print('  residual R1 vanishes (gadget implements doubling): %d / %d'%(ok,tested))
print('  residual R1 nonzero  (chord only, doubling NOT representable): %d / %d'%(nz,tested))
