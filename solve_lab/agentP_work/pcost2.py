#!/usr/bin/env python3
"""Agent P: corrected costing. For each law-block, its 3 law-congruence atoms live
contiguously in the congruence section (pos>19000) together with their handle-defs.
Cost a corruption of that block's law output = union of equations touched."""
import pickle,sys
from collections import Counter,defaultdict
sys.set_int_max_str_digits(10**7)
W='/home/user/integer_solver/solve_lab/agentP_work/'
D=pickle.load(open(W+'model4.pkl','rb')); AP=D['AP']; rows=D['rows']
S=pickle.load(open(W+'slp.pkl','rb')); topo=S['topo']; outof=S['outof']
B=pickle.load(open(W+'blocks.pkl','rb'))
T=pickle.load(open(W+'topo.pkl','rb')); supp=T['supp']; src=T['src']
pos={a:i for i,a in enumerate(topo)}
cover=defaultdict(set)
for ei,r in enumerate(rows):
    for c,a in r['row']: cover[a].add(ei)
v2a=defaultdict(list)
for i,ap in enumerate(AP):
    for x in set(y for m in ap for y in m): v2a[x].append(i)

blockcong={}
for j,b in enumerate(B):
    trip=[]
    for (ca,cb,tc) in b['outs']:
        cong=[a for a in v2a[tc] if pos[a]>19000]
        pack=[]
        for a in cong:
            hd=[]
            for x in set(y for m in AP[a] for y in m):
                for h in v2a[x]:
                    if outof[h]==x and h!=a and pos[h]>19000: hd.append(h)
            pack.append((a,hd))
        trip.append(pack)
    blockcong[j]=trip
bad=[j for j in range(len(B)) if any(len(t)!=1 for t in blockcong[j])]
print("blocks whose 3 O_k each map to exactly 1 congruence atom:",len(B)-len(bad),"of",len(B))

rank=[]
for j in range(len(B)):
    if j in bad: continue
    opts=[]
    for w in ((0,1),(0,2),(1,2)):
        for mode in ('cong_only','cong+handle'):
            A=set()
            for k in w:
                a,hd=blockcong[j][k][0]
                A.add(a)
                if mode=='cong+handle': A.update(hd)
            eq=set()
            for a in A: eq|=cover[a]
            opts.append((len(eq),len(A),mode,w,sorted(pos[a] for a in A)))
    opts.sort()
    rank.append((opts[0][0],opts[0][1],j,opts[0][2],opts[0][3],opts[0][4]))
rank.sort()
print()
print("cheapest placements (touched equations / atoms broken):")
for r in rank[:18]:
    j=r[2]
    kinds=tuple(k[0] for k in src[j])
    print("  touched=%2d atoms=%d block=%3d %-11s break=%s supp=%3d inputs=%s pos=%s"%(
        r[0],r[1],j,r[3],r[4],len(supp[j]),kinds,r[5]))
print()
# where does the deliverable's own block sit?
dl=[a for a in range(len(AP)) if pos[a] in (36292,36294)]
print("deliverable law-congruence atoms at slp 36292,36294 ->",dl)
for j in range(len(B)):
    if j in bad: continue
    for k in range(3):
        if pos[blockcong[j][k][0][0]] in (36292,36294):
            print("   they belong to block",j,"supp",len(supp[j]),"inputs",tuple(x[0] for x in src[j]))
pickle.dump({'blockcong':blockcong,'rank':rank},open(W+'cost2.pkl','wb'))
