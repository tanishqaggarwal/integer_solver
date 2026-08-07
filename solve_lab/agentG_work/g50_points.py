"""Pair up each message bit's load constants into secp256k1 POINTS, then test whether
the 256 points form a doubling chain (=> the message is a scalar and the instance is
a scalar multiplication / ECDLP) or an unstructured table (=> enumerable)."""
import os, sys, pickle, collections, itertools
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentG_work')
import gpt
from gsym2 import L, ad, P
d=pickle.load(open('/home/user/integer_solver/solve_lab/agentG_work/loads.pkl','rb'))
tab=d['tab']
n=gpt.n
def tox(H): return ((H+gpt.K*pow(3,-1,P))%P)*gpt.IU2%P
def toy(H): return H*gpt.IU3%P
pts={}
for b,ls in tab.items():
    cands=set()
    consts=[(xB,H) for _,xB,H in ls]
    for (xa,Ha),(xb2,Hb) in itertools.permutations(consts,2):
        X=tox(Ha); Y=toy(Hb)
        if (Y*Y-pow(X,3,P)-7)%P==0: cands.add((X,Y))
    if cands: pts[b]=sorted(cands)
print('message bits with a well-formed curve point among their load constants: %d of %d'%(len(pts),len(tab)))
mult=[b for b in pts if len(pts[b])>1]
print('bits with more than one point: %d'%len(mult))
bs=sorted(pts)
for b in bs[:8]: print('  x%-6d %d point(s): %s'%(b,len(pts[b]),str(pts[b][0])[:70]))
# doubling-chain test
print('\n--- doubling / small-multiple structure ---')
allp=[(b,q) for b in bs for q in pts[b]]
print('total distinct points: %d'%len(set(q for _,q in allp)))
idx={q:b for b,q in allp}
def dbl(Q): return gpt.add(Q,Q)
hits=0
for b,q in allp:
    D=dbl(q)
    if D in idx: hits+=1; print('   [2]P(x%d) = P(x%d)'%(b,idx[D]))
print('doubling hits: %d'%hits)
# is any point a small multiple of G?
Gp=(gpt.Gx,gpt.Gy); R=None; small={}
for k in range(1,20001):
    R=gpt.add(R,Gp); small[R]=k
sm=[(b,small[q]) for b,q in allp if q in small]
print('points that are [k]G for k<=20000:',sm[:20])
# differences between consecutive bits
print('\n--- pairwise relation of the first few points ---')
sel=[allp[i] for i in range(min(6,len(allp)))]
for (b1,q1),(b2,q2) in itertools.combinations(sel,2):
    dd=gpt.sub(q2,q1)
    print('   P(x%d)-P(x%d) = %s'%(b2,b1,str(dd)[:50]))
pickle.dump(pts,open('/home/user/integer_solver/solve_lab/agentG_work/bitpoints.pkl','wb'))
