#!/usr/bin/env python3
"""Q-11e: the liveness graph.  Each slot has live bits (s1,s2) and emits OR(s1,s2).  If every slot's
live bits are its children's OR outputs (or leaf selectors, or hard zeros), the slots form a tree."""
import pickle,collections,re,json
exec(open('qmux2.py').read().split('ST=[x for x')[0])
ST=[x for x in json.load(open('qstages.json'))['stages'] if 'u3' in x]
leaf={int(g):v for g,v in json.load(open('qleaf.json')).items()}
LEAFSEL={int(g) for g in leaf}
SUBP=re.compile(r'^x_(\d+) - \(x_(\d+) - x_(\d+)\)$')
ZERO=re.compile(r'^x_(\d+) - 0$')
ORof={}; zeros=set()
for s,vs in terms:
    m=SUBP.match(s)
    if m:
        t,u,v=map(int,m.groups())
        if u in summ and v in prod and set(summ[u])==set(prod[v]): ORof[t]=tuple(sorted(map(root,summ[u])))
        continue
    m=ZERO.match(s)
    if m: zeros.add(int(m.group(1)))
sl={}
for gi,g in enumerate(ST):
    for P,Q in (((g['ua'],g['ub'],g['u3']),(g['ya'],g['yb'],g['y3'])),
                ((g['ya'],g['yb'],g['y3']),(g['ua'],g['ub'],g['u3'])),
                ((g['ua'],g['ub'],g['y3']),(g['ya'],g['yb'],g['u3'])),
                ((g['ya'],g['yb'],g['u3']),(g['ua'],g['ub'],g['y3']))):
        got=None
        for cA,cB,cC,X in mux(*P):
            q=quad(cA,cB,cC)
            if not q: continue
            for a2,b2,c2,Y in mux(*Q):
                if (cA,cB,cC)==(a2,b2,c2): got=q; break
            if got: break
        if got: break
    sl[gi]=tuple(map(root,got))
# the OR that this slot emits
emit={}
for t,(a,b) in ORof.items():
    for gi,(s1,s2) in sl.items():
        if tuple(sorted((s1,s2)))==tuple(sorted((a,b))): emit.setdefault(gi,root(t))
print('slots %d ; slots emitting an OR of their own two live bits: %d'%(len(ST),len(emit)))
child_of={}
for gi,(s1,s2) in sl.items():
    for s in (s1,s2):
        for gj,e in emit.items():
            if e==s: child_of.setdefault((gi,s),gj)
edges=[(gj,gi) for (gi,s),gj in child_of.items()]
print('parent<-child edges recovered: %d'%len(edges))
par=collections.defaultdict(list); haspar=set()
for c,p in edges: par[p].append(c); haspar.add(c)
roots=[gi for gi in sl if gi not in haspar]
print('slots with no parent (tree roots): %d'%len(roots))
kinds=collections.Counter()
for gi,(s1,s2) in sl.items():
    for s in (s1,s2):
        kinds['leaf selector' if s in LEAFSEL else ('hard zero' if s in zeros else ('child OR' if any(e==s for e in emit.values()) else 'other'))]+=1
print('live-bit sources:',dict(kinds))
# reachability from the root
if len(roots)==1:
    seen=set(); stk=[roots[0]]
    while stk:
        x=stk.pop()
        if x in seen: continue
        seen.add(x); stk+=par.get(x,[])
    print('slots reachable from the single root: %d / %d'%(len(seen),len(ST)))
    lf=set()
    for gi in seen:
        for s in sl[gi]:
            if s in LEAFSEL: lf.add(s)
    print('distinct leaf selectors under the root: %d / 256'%len(lf))
