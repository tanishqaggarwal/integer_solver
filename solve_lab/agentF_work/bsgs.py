#!/usr/bin/env python3
import json,sys,time
d=json.load(open('curve.json'))
p=115792089237316195423570985008687907853269984665640564039457584007908834671663
b=int(d['b']); N=int(d['N'])
T=(int(d['T'][0]),int(d['T'][1]))
P0=tuple(int(z) for z in d['pts'][str(d['roots'][0])])
def add(P,Q):
    if P is None: return Q
    if Q is None: return P
    x1,y1=P;x2,y2=Q
    if x1==x2:
        if (y1+y2)%p==0: return None
        l=3*x1*x1%p*pow(2*y1,p-2,p)%p
    else: l=(y2-y1)*pow(x2-x1,p-2,p)%p
    x3=(l*l-x1-x2)%p; return (x3,(l*(x1-x3)-y1)%p)
def mul(k,P):
    k%=N; R=None; Q=P
    while k:
        if k&1: R=add(R,Q)
        Q=add(Q,Q); k>>=1
    return R
def neg(P): return None if P is None else (P[0],(-P[1])%p)
# quick structured checks
print('quick structured-log checks',flush=True)
for k in list(range(1,5000))+[2**i for i in range(1,256)]+[N-1,(N-1)//2,(N+1)//2,3,7]:
    if mul(k,P0)==T: print('FOUND small/structured log k=',k,flush=True); sys.exit()
print('  none up to 5000 / powers of two',flush=True)
M=1<<20
t0=time.time()
baby={}
Q=None
for j in range(M):
    key = 0 if Q is None else Q[0]&((1<<62)-1)
    baby.setdefault(key,[]).append(j)
    Q=add(Q,P0)
    if j%200000==0: print('  baby',j,round(time.time()-t0,1),flush=True)
print('baby done',M,round(time.time()-t0,1),flush=True)
S=mul(M,P0); Sn=neg(S)
cur=T
for i in range(1<<20):
    key = 0 if cur is None else cur[0]&((1<<62)-1)
    if key in baby:
        for j in baby[key]:
            if mul(j,P0)==cur:
                k=(i*M+j)%N
                print('FOUND LOG k =',k,flush=True); sys.exit()
    cur=add(cur,Sn)
    if i%100000==0: print('  giant',i,round(time.time()-t0,1),flush=True)
print('NO LOG below 2^40',flush=True)
