#!/usr/bin/env python3
"""Agent X: rebuild + independently verify the curve, the target T and the ladder 2^i*G,
straight from agent Q's raw instance-derived files (read-only).  Writes xdata.json."""
import json, os, sys
HERE=os.path.dirname(os.path.abspath(__file__))
Q=os.path.join(HERE,'..','agentQ_work')
p = 115792089237316195423570985008687907853269984665640564039457584007908834671663
K = 97553848499418123410591666447050222001188385549510401465815187079080512838891
N = 115792089237316195423570985008687907852837564279074904382605163141518161494337
cs = K*pow(3,p-2,p)%p
cur=json.load(open(os.path.join(Q,'curve.json')))
A_=int(cur['a']); B_=int(cur['b'])
assert int(cur['p'])==p and int(cur['c_shift'])==cs
print('p == 2^256-2^32-977 (secp256k1 prime):', p==2**256-2**32-977)
print('N == secp256k1 group order:', N==0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141)
def inv(z): return pow(z,p-2,p)
def add(P,Q_):
    if P is None: return Q_
    if Q_ is None: return P
    x1,y1=P; x2,y2=Q_
    if (x1-x2)%p==0:
        if (y1+y2)%p==0: return None
        l=(3*x1*x1+A_)%p*inv(2*y1%p)%p
    else:
        l=(y2-y1)%p*inv((x2-x1)%p)%p
    x3=(l*l-x1-x2)%p
    return (x3,(l*(x1-x3)-y1)%p)
def neg(P): return None if P is None else (P[0],(-p+p-P[1])%p if False else (-P[1])%p)
def mul(k,P):
    if k<0: k,P=-k,neg(P)
    R=None;Qq=P
    while k>0:
        if k&1: R=add(R,Qq)
        Qq=add(Qq,Qq); k>>=1
    return R
def oncur(P): return P is None or (P[1]*P[1]-pow(P[0],3,p)-A_*P[0]-B_)%p==0

leaf={int(g):(int(v[0]),int(v[1])) for g,v in json.load(open(os.path.join(Q,'qleaf.json'))).items()}
print('leaves: %d, all on curve: %s'%(len(leaf), all(oncur(v) for v in leaf.values())))
lad=json.load(open(os.path.join(Q,'qladder.json')))
e2s={int(k):int(v) for k,v in lad['exp2sel'].items()}
print('ladder exponents: %d (0..%d)'%(len(e2s),max(e2s)))
G=leaf[e2s[0]]
# independent check: L_i == 2^i*G by repeated doubling
cur_pt=G; ok=0; bad=[]
for i in range(256):
    if i in e2s:
        if leaf[e2s[i]]==cur_pt: ok+=1
        else: bad.append(i)
    cur_pt=add(cur_pt,cur_pt)
print('ladder L_i == 2^i*G :  %d/%d  bad=%s'%(ok,len(e2s),bad))
assert mul(N,G) is None, 'N*G != O'
print('N*G == O : True   (G has prime order N)')
# target
C1=91416258160755509149180373473728639746431157665678710450404458852172057265575180278101002
C2=125787314747601108116039725163361763116550465675981151838811516827327919228823597744635626
T=((C1+cs)%p,C2%p)
print('T on curve:', oncur(T))
assert mul(N,T) is None
print('N*T == O : True')
# cross-check against Q's hardcoded values in lowwt.py
assert G==(31917591553801470078828036568057743875467637605644620066197178005619323650152,
           83364444556352143115103874010002344754157095926378075484791050960431190202517), 'G mismatch'
assert T==(30121525689829097248416773597728729849687459852468451992398421980273013515302,
           44859544763832475231923253825569092119321525945631045653619508440821028887), 'T mismatch'
print('G and T match the values Q searched with: True')
# full ladder as points, index i -> 2^i*G
P=[G]
for i in range(255): P.append(add(P[-1],P[-1]))
assert all(oncur(q) for q in P)
json.dump({'p':str(p),'a':str(A_),'b':str(B_),'N':str(N),
           'G':[str(G[0]),str(G[1])],'T':[str(T[0]),str(T[1])],
           'ladder':[[str(x),str(y)] for x,y in P],
           'exp2sel':{str(k):v for k,v in e2s.items()}},
          open(os.path.join(HERE,'xdata.json'),'w'))
print('wrote xdata.json  (256 ladder points)')
