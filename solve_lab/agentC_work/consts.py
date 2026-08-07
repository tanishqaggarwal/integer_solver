import sys, re, json, time
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentC_work')
from ec import *
chain=[int(b) for b in json.load(open('/home/user/integer_solver/solve_lab/agentC_work/chain.json'))]
pts=leafpoints()
G=pts[chain[0]]
C=json.load(open('/home/user/integer_solver/solve_lab/agentC_work/curve.json'))
Q=(int(C['Q'][0]),int(C['Q'][1]))
print('chain[255]==2^255G?', mul(pow(2,255),G)==pts[chain[255]])
print('sanity: G*1',mul(1,G)==G)
# all integer literals in the file
t=time.time()
S=set()
with open('/home/user/integer_solver/EQUATIONS.txt') as f:
    for line in f:
        for m in re.finditer(r'(?<![x_\d])(\d{6,})', line):
            S.add(int(m.group(1)))
print('distinct literals >=6 digits:',len(S),time.time()-t)
big=[c for c in S if c>10**20]
print('literals > 1e20:',len(big))
cands=set()
for c in big:
    cands.add(c); cands.add(c%N)
    cands.add((-c)%N)
P2_=2**256-2**32-977
cands.add(P2_); cands.add(N)
print('scalar candidates',len(cands))
t=time.time(); hit=None
for i,k in enumerate(cands):
    if mul(k,G)==Q: hit=k; print('*** DLOG FOUND (literal):',k); break
    if i%200==0: print('  ',i,time.time()-t,flush=True)
print('done',hit,time.time()-t)
json.dump({'hit':str(hit)},open('/home/user/integer_solver/solve_lab/agentC_work/consts_hit.json','w'))
