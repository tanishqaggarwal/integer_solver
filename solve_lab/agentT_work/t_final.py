#!/usr/bin/env python3
"""AUDIT T22 -- (a) are the 3 omitted-but-incident atoms genuine handles?  (b) decode a 278-path."""
import os,sys,pickle,collections,re
T=os.path.dirname(os.path.abspath(__file__)); LAB=os.path.abspath(os.path.join(T,'..'))
F=os.path.join(LAB,'agentF_work'); sys.path.insert(0,F); sys.path.insert(0,LAB)
sys.path.insert(0,os.path.join(LAB,'agentE_work'))
from circ2 import vars_of
import checker as CK, engine as E
d=pickle.load(open(os.path.join(F,'circ4.pkl'),'rb'))
atoms=d['atoms']; eqrows=d['eqrows']; names=list(atoms); idx={a:i for i,a in enumerate(names)}
avars=[frozenset(vars_of(atoms[a])) for a in names]
v2a=collections.defaultdict(set)
for i,vs in enumerate(avars):
    for u in vs: v2a[u].add(i)
a2e=collections.defaultdict(set)
for e,row in enumerate(eqrows):
    for k,a in row: a2e[idx[a]].add(e)
codes,varsets=CK.load_equations()
v2e=collections.defaultdict(set)
for e,vs in enumerate(varsets):
    for u in vs: v2e[u].add(e)
M25=set([5324,9041,11226,15558,21000,22534,22997,28929,29330,32026,35512,38051,
     12231,12270,12350,14584,18673,22044,29125])
print('== (a) the 3 atoms L omits that ARE incident to the baseline-failing set ==')
for u in [10422,15120,35531]:
    A=v2a[u]
    ok = len(A)==1
    a=next(iter(A))
    match = v2e[u]==a2e[a]
    print('  u=x%-6d free=%-5s  in %d atom(s)  eqs(u)==eqs(atom): %-5s'%(u,E.definer[u] is None,len(A),match))
    print('      atom : %s'%names[a][:78])
    h=re.match(r'^\(x(\d+)-',names[a].replace(' ',''))
    hv=int(h.group(1)) if h else None
    g=[names[j] for j in v2a[hv] if j!=a] if hv else []
    print('      guard: %s'%(g[0][:78] if g else '(none)'))
    print('      touches baseline-failing equations: %s'%sorted(v2e[u]&M25))
print('\n  => these satisfy the SAME criterion L verified on its 3,681 (free, unique atom, eqs match).')
print('\n== (b) a concrete multi-hop alias path ==')
Mo=pickle.load(open(os.path.join(LAB,'agentL_work','full_model.pkl'),'rb'))
OUT=Mo['OUT']; tree=Mo['tree']
pairs=[]
for n in tree:
    if n not in OUT or not tree[n] or len(tree[n])!=2: continue
    ca,cb=tree[n]
    for j,slot in enumerate(OUT[n]):
        for side,ch in (('va',ca),('vb',cb)):
            w=slot.get(side)
            if ch not in OUT or w is None: continue
            co=OUT[ch][j].get('out')
            if co is None or w==co: continue
            pairs.append((n,side,ch,w,co))
multi=[t for t in pairs if not (v2a[t[3]]&v2a[t[4]])]
def path(w,co,maxd=4):
    seen={w}; fr=[(w,[])]
    for dep in range(maxd):
        nf=[]
        for x,pa in fr:
            for i in v2a[x]:
                for y in avars[i]:
                    if y in seen: continue
                    if y==co: return pa+[i]
                    seen.add(y); nf.append((y,pa+[i]))
        fr=nf
    return None
for n,side,ch,w,co in multi[:3]:
    p=path(w,co)
    print('\n  node x%d.%s  <- child x%d   (parent slot x%d, child out x%d)'%(n,side,ch,w,co))
    for i in p: print('      %s'%names[i][:88])
