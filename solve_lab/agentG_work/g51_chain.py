"""Order the 256 message-bit points into the doubling chain and identify its root."""
import os, sys, pickle
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentG_work')
import gpt
from gsym2 import L, ad, P
pts=pickle.load(open('/home/user/integer_solver/solve_lab/agentG_work/bitpoints.pkl','rb'))
allp=[(b,q) for b in sorted(pts) for q in pts[b]]
idx={q:b for b,q in allp}
succ={}; pred={}
for b,q in allp:
    D=gpt.add(q,q)
    if D in idx: succ[b]=idx[D]; pred[idx[D]]=b
roots=[b for b,_ in allp if b not in pred]
print('chain roots (no predecessor): %s'%roots)
print('chain leaves (no successor) : %s'%[b for b,_ in allp if b not in succ])
order=[]; cur=roots[0] if roots else allp[0][0]
seen=set()
while cur is not None and cur not in seen:
    seen.add(cur); order.append(cur); cur=succ.get(cur)
print('chain length: %d of %d points'%(len(order),len(allp)))
P0=pts[order[0]][0]
print('\nroot point P0 = %s'%(P0,))
Gp=(gpt.Gx,gpt.Gy)
print('P0 == G ?', P0==Gp)
print('P0 == -G ?', P0==(Gp[0],(-Gp[1])%P))
print('[n]P0 = %s'%('O' if gpt.mul(gpt.n,P0) is None else 'NOT O'))
# verify the chain really is [2^i]P0
bad=0
Q=P0
for i,b in enumerate(order):
    if pts[b][0]!=Q: bad+=1
    Q=gpt.add(Q,Q)
print('positions where P(bit_i) != [2^i]P0 : %d of %d'%(bad,len(order)))
# is P0 a small multiple of G?
R=None; small={}
for k in range(1,200001):
    R=gpt.add(R,Gp); small[R]=k
print('P0 = [k]G for k<=200000 ?', small.get(P0))
# relation of the base points P1,P2,P3 to the chain
import g46_table as T
base=T.frame([])
B1,B2,B3=base['pts']
for nm,Q in [('P1',B1),('P2',B2),('P3',B3)]:
    print('%s in chain? %s ; %s = [k]G small? %s'%(nm, idx.get(Q), nm, small.get(Q)))
D=base['D']
print('D = P3-(P1+P2) in chain? %s ; small multiple of G? %s'%(idx.get(D),small.get(D)))
pickle.dump({'order':order,'P0':P0},open('/home/user/integer_solver/solve_lab/agentG_work/chain.pkl','wb'))
