#!/usr/bin/env python3
"""Greedy topological scheduling of definition atoms."""
import sys,os,pickle,collections,time
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE)
from circ2 import vars_of
d=pickle.load(open(os.path.join(HERE,'circ4.pkl'),'rb'))
atoms,cls,eqrows=d['atoms'],d['cls'],d['eqrows']
NV=38748
defs=collections.defaultdict(list); cons=[]
for s,c in cls.items():
    if c[0]=='def': defs[c[1]].append(s)
    else: cons.append(s)
# dependency
rhsvars={}
for s,c in cls.items():
    if c[0]=='def': rhsvars[s]=vars_of(c[2])
    else: rhsvars[s]=vars_of(c[2])
known=set(range(NV))-set(defs)          # free vars
print('free',len(known))
pend=collections.defaultdict(list)      # var -> list of atoms waiting on it
ready=[]
cnt={}
for v,ss in defs.items():
    for s in ss:
        miss=[u for u in rhsvars[s] if u not in known]
        cnt[s]=len(miss)
        if not miss: ready.append(s)
        else:
            for u in set(miss): pend[u].append(s)
order=[]; usedcons=[]
while ready:
    s=ready.pop()
    v=cls[s][1]
    if v in known:
        usedcons.append(s); continue
    known.add(v); order.append(s)
    for t in pend.get(v,[]):
        cnt[t]-=1
        if cnt[t]==0: ready.append(t)
print('scheduled defs',len(order),'known vars',len(known),'redundant-def-as-constraint',len(usedcons))
blocked=[s for s in cnt if cnt[s]>0 and cls[s][1] not in known]
print('blocked (cyclic) defs',len([s for s in cnt if cnt[s]>0]))
unk=set(range(NV))-known
print('unknown vars after schedule',len(unk))
pickle.dump({'order':order,'usedcons':usedcons,'cons':cons,'unk':sorted(unk)},open(os.path.join(HERE,'sched.pkl'),'wb'))
