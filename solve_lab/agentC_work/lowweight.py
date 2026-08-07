"""Meet in the middle for low Hamming weight k: subsets of size <= 6."""
import sys, json, time, itertools
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentC_work')
from ec import *
chain=[int(b) for b in json.load(open('/home/user/integer_solver/solve_lab/agentC_work/chain.json'))]
pts=leafpoints()
G=pts[chain[0]]
C=json.load(open('/home/user/integer_solver/solve_lab/agentC_work/curve.json'))
Q=(int(C['Q'][0]),int(C['Q'][1]))
Pl=[pts[b] for b in chain]
n=len(Pl)
T0=time.time()
# size 1
for i,pp in enumerate(Pl):
    if pp==Q: print('WEIGHT1 SOLUTION bit',i)
# pairs
pairsum={}
for i in range(n):
    for j in range(i+1,n):
        s=add(Pl[i],Pl[j])
        if s is not None: pairsum.setdefault(s,(i,j))
print('pair sums',len(pairsum),time.time()-T0,flush=True)
if Q in pairsum: print('WEIGHT2 SOLUTION',pairsum[Q])
# weight 3: Q - Pl[k] in pairsum
for k in range(n):
    r=add(Q,neg(Pl[k]))
    if r in pairsum: print('WEIGHT3 SOLUTION',k,pairsum[r])
print('w3 done',time.time()-T0,flush=True)
# weight 4: Q - pairsum in pairsum
for s,(i,j) in pairsum.items():
    r=add(Q,neg(s))
    if r in pairsum: print('WEIGHT4 SOLUTION',(i,j),pairsum[r])
print('w4 done',time.time()-T0,flush=True)
# triples
tri={}
cnt=0
for i in range(n):
    for j in range(i+1,n):
        s2=add(Pl[i],Pl[j])
        for k in range(j+1,n):
            s=add(s2,Pl[k])
            if s is not None: tri[s]=(i,j,k)
        cnt+=1
    if i%32==0: print('  tri',i,len(tri),time.time()-T0,flush=True)
print('triple sums',len(tri),time.time()-T0,flush=True)
for s,t in tri.items():
    r=add(Q,neg(s))
    if r in pairsum: print('WEIGHT5 SOLUTION',t,pairsum[r])
    if r in tri: print('WEIGHT6 SOLUTION',t,tri[r])
print('done',time.time()-T0)
