import sys, json
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentC_work')
from ec import *
from ort import *
pts=leafpoints()
bits=sorted(pts)
S=set(pts.values())
print('distinct leaf points:',len(S),'of',len(pts))
dbl_hits=sum(1 for b in bits if add(pts[b],pts[b]) in S)
print('leaf points whose DOUBLE is also a leaf point:',dbl_hits,'/',len(bits))
# side membership
sides={}
for r,tag in [(8599,'s1a'),(21839,'s1b'),(25956,'s2a'),(7304,'s2b')]:
    for x in leaves(r):
        if x in pts: sides[x]=tag
import collections
print(collections.Counter(sides.values()))
# build chain
rev={v:k for k,v in pts.items()}
succ={}
for b in bits:
    d=add(pts[b],pts[b])
    if d in rev: succ[b]=rev[d]
print('chain edges',len(succ))
# find chain roots
tgt=set(succ.values())
roots=[b for b in bits if b not in tgt]
print('chain roots',len(roots))
if len(succ)>200:
    r=roots[0] if roots else bits[0]
    ch=[r]; seen={r}
    while ch[-1] in succ and succ[ch[-1]] not in seen:
        ch.append(succ[ch[-1]]); seen.add(ch[-1])
    print('longest chain from root length',len(ch))
    json.dump([str(b) for b in ch],open('/home/user/integer_solver/solve_lab/agentC_work/chain.json','w'))
    print('chain sides:',[sides.get(b) for b in ch[:20]])
