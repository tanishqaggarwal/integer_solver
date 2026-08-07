#!/usr/bin/env python3
"""Agent P: cost (in touched equations) of corrupting each stage's law output,
and of corrupting each leaf's coordinate load."""
import pickle,sys,json
from collections import Counter,defaultdict
sys.set_int_max_str_digits(10**7)
W='/home/user/integer_solver/solve_lab/agentP_work/'
D=pickle.load(open(W+'model4.pkl','rb')); AP=D['AP']; rows=D['rows']
S=pickle.load(open(W+'slp.pkl','rb')); topo=S['topo']; outof=S['outof']
B=pickle.load(open(W+'blocks.pkl','rb'))
LEAVES=pickle.load(open(W+'leaves.pkl','rb'))
pos={a:i for i,a in enumerate(topo)}
cover=defaultdict(set)
for ei,r in enumerate(rows):
    for c,a in r['row']: cover[a].add(ei)

# var -> atoms containing it
v2a=defaultdict(list)
for i,ap in enumerate(AP):
    for x in set(y for m in ap for y in m): v2a[x].append(i)

def congr_for(var):
    """atoms other than the SLP-definition of var that contain var"""
    dp=None
    out=[]
    for a in v2a[var]:
        if outof[a]==var and pos[a]<31000: dp=a
        else: out.append(a)
    return dp,out

recs=[]
for j,b in enumerate(B):
    grp=[]
    for (ca,cb,tc) in b['outs']:
        dp,rest=congr_for(tc)
        # each rest atom is  O*ab - c*handle ; find the handle's own def atom
        pack=[]
        for a in rest:
            hs=[x for m in AP[a] for x in m if x!=tc]
            hdef=[h for x in hs for h in v2a[x] if outof[h]==x and h!=a]
            pack.append((a,hdef))
        grp.append((tc,pack))
    recs.append(grp)

def atomset(j, which):
    """atoms broken if we break congruences `which` (indices into outs) of block j,
       taking both the congruence atom and its handle-def atom."""
    S=set()
    for k in which:
        tc,pack=recs[j][k]
        for a,hd in pack:
            S.add(a); S.update(hd)
    return S

print("cost of corrupting one stage's law output (break 2 of its 3 congruences):")
res=[]
for j in range(len(B)):
    best=None
    for w in ((0,1),(0,2),(1,2)):
        S=atomset(j,w)
        eq=set()
        for a in S: eq|=cover[a]
        if best is None or len(eq)<best[0]: best=(len(eq),w,len(S),sorted(pos[a] for a in S))
    res.append((best[0],best[2],j,best[1],best[3]))
res.sort()
for r in res[:15]:
    print("  touched=%3d  atoms=%d  block=%3d  break=%s  slp_pos=%s"%r)
print("  ... worst:",res[-1][:4])
pickle.dump({'recs':recs,'res':res},open(W+'cost.pkl','wb'))

print()
print("cost of corrupting one leaf's two coordinate loads:")
lres=[]
byleaf=defaultdict(list)
for (a,bb,k) in LEAVES: byleaf[a].append(bb)
for selv,coords in byleaf.items():
    S=set()
    for c in coords:
        for a in v2a[c]:
            if pos[a]>=19000:      # the congruence section
                S.add(a)
                for x in set(y for m in AP[a] for y in m):
                    for h in v2a[x]:
                        if outof[h]==x and h!=a and pos[h]>=19000: S.add(h)
    eq=set()
    for a in S: eq|=cover[a]
    lres.append((len(eq),len(S),selv,sorted(pos[a] for a in S)))
lres.sort()
for r in lres[:10]: print("  touched=%3d atoms=%2d selector=x%d slp_pos=%s"%r)
print("  ... worst:",lres[-1][:3])
