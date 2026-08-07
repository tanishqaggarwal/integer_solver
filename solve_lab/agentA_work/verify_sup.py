import sys, json, random; sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentA_work')
from fractions import Fraction as F
import env, lib as L
from agrow import model
import amk_model as MK
A,K,R,rows,QUAD=model([37887,41906])
aff=[(e,c,lin) for e,c,lin,hq in rows if not hq]
nk=len(K)
NZ=[(e,c,lin) for e,c,lin in aff if lin]
n=len(NZ); EQ=[e for e,_,_ in NZ]
N=[[lin.get(j,0) for j in range(nk)] for e,c,lin in NZ]
def rankQ(idx):
    M=[[F(N[i][j]) for j in range(nk)] for i in idx]
    r=0
    for c in range(nk):
        pr=None
        for i in range(r,len(M)):
            if M[i][c]!=0: pr=i;break
        if pr is None: continue
        M[r],M[pr]=M[pr],M[r]
        pv=M[r][c]; M[r]=[x/pv for x in M[r]]
        for i in range(len(M)):
            if i!=r and M[i][c]!=0:
                f=M[i][c]; M[i]=[x-f*y for x,y in zip(M[i],M[r])]
        r+=1
        if r==len(M): break
    return r
print('rank(all %d rows) over Q = %d'%(n,rankQ(range(n))))
sups=[set(s) for s in json.load(open('/home/user/integer_solver/solve_lab/agentA_work/supports.json'))]
random.seed(0)
ok=0; bad=0
for s in random.sample(sups,min(25,len(sups))):
    Z=[i for i in range(n) if EQ[i] not in s]
    r=rankQ(Z)
    if r<nk: ok+=1
    else: bad+=1; 
    print('  D=%s rank(complement)=%d %s'%(sorted(s),r,'OK' if r<nk else 'FALSE POSITIVE'))
print('verified %d real, %d false positives'%(ok,bad))
