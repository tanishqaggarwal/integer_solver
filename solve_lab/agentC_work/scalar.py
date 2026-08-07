import sys, json, collections
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentC_work')
from ec import *
from ort import *
pts=leafpoints()
chain=[int(b) for b in json.load(open('/home/user/integer_solver/solve_lab/agentC_work/chain.json'))]
idx={b:i for i,b in enumerate(chain)}
bs=set(chain)
# every atom containing >=4 leaf bits
cnt=collections.Counter()
for a in range(L.NA):
    k=len(L.avars[a]&bs)
    if k>=4: cnt[a]=k
print('atoms with >=4 leaf bits:',len(cnt))
for a,k in cnt.most_common(15):
    print('a%-6d nbits=%-3d out=%s'%(a,k,L.atom_out.get(a)))
    print('   ',L.atom_src[a][:300])
# also: check atoms where a leaf bit has coefficient that is a power of two
print()
pw={1<<i:i for i in range(300)}
hits=[]
for a in range(L.NA):
    Pp=L.polys[a]
    tw=[(m[0],c) for m,c in Pp.items() if len(m)==1 and m[0] in bs and abs(c) in pw]
    if len(tw)>=3: hits.append((a,tw))
print('atoms with >=3 leaf bits at power-of-two coefficients:',len(hits))
for a,tw in hits[:5]:
    print('a%d'%a,[(u,c) for u,c in tw][:8]); print('   ',L.atom_src[a][:250])
