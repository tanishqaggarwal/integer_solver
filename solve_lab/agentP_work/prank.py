#!/usr/bin/env python3
"""Agent P: the 6-parameter rank computation at block 2, |S|=2.

Parameters (all shifts are multiples of P, written  i_k -> i_k + P*mu_k):
  i1..i4  leaf coordinates, legal step c_leaf*P  =>  mu_k ranges over c_leaf*Z
  i5,i6   block law outputs, legal step P        =>  mu_k ranges over Z
Exact expansion (sN1=+1, sN2=-1 as extracted):
  N1 = E*A^2 - B^2 ,  N2 = A*(i3+i6) - B*(i2-i5)
  A = i1-i2 , B = i4-i3 , E = i1+i2+i5+Q , H = i3+i6 , J = i2-i5
  a = mu1+mu2+mu5 , b = mu1-mu2 , d = mu4-mu3 , g = mu3+mu6 , h2 = mu2-mu5
  n1' = n1 + (2*E*A*b + a*A^2 - 2*B*d) + P*(E*b^2 + 2*a*A*b - d^2) + P^2*a*b^2
  n2' = n2 + (A*g + H*b - A*h2 - J*d) + P*(b*g - d*h2)
Condition per k:  c_k | c_k1*n1' + c_k2*n2'
"""
import pickle,sys
from math import gcd
from itertools import combinations
sys.set_int_max_str_digits(10**7)
W='/home/user/integer_solver/solve_lab/agentP_work/'
D=pickle.load(open(W+'model4.pkl','rb')); AP=D['AP']
S=pickle.load(open(W+'slp.pkl','rb')); topo=S['topo']; outof=S['outof']
B=pickle.load(open(W+'blocks.pkl','rb'))
import plift5, pfold as F
P=plift5.P; Q=F.Q
pos={a:i for i,a in enumerate(topo)}

def fac(n):
    f={};d=2
    while d*d<=n:
        while n%d==0: f[d]=f.get(d,0)+1; n//=d
        d+=1
    if n>1: f[n]=f.get(n,0)+1
    return f

# configuration: the first leaf-leaf merge
for j,row in enumerate(F.SRC):
    if all(k[0]=='L' for k in row): sel={row[0][1],row[1][1]}; jj=j; break
val,obs,und,sz,live=plift5.build(sel)
b=B[jj]
i1,i2,i3,i4,i5,i6=[val[b[k]] for k in ('i1','i2','i3','i4','i5','i6')]
A=i1-i2; Bv=i4-i3; E=i1+i2+i5+Q; H=i3+i6; J=i2-i5
N1=E*A*A-Bv*Bv; N2=A*H-Bv*J
assert N1%P==0 and N2%P==0, "base point is not mod-P valid"
n1=N1//P; n2=N2//P
print("block %d, |S|=2, base point verified: N1 and N2 both divisible by P"%jj)

# legal step for each of i1..i4 (leaf coordinates): x = K + c*P*u  -> mu step = c
step=[1]*6
for idx,vn in enumerate(('i1','i2','i3','i4')):
    v=b[vn]
    cs=[]
    for a in range(len(AP)):
        if pos[a]<19000: continue
        ap=AP[a]
        vs={y for m in ap for y in m}
        if v in vs and any(x in plift5.HANDLE for x in vs):
            for m,c in ap.items():
                if len(m)==1 and m[0] in plift5.HANDLE: cs.append(abs(c))
    step[idx]=cs[0] if cs else 1
print("legal mu-steps for (i1,i2,i3,i4,i5,i6):",step)

rows3=[]
for k,(ca,cb,tc) in enumerate(b['outs']):
    cong=[a for a in range(len(AP)) if tc in {y for m in AP[a] for y in m} and pos[a]>19000]
    ck=1
    for m,c in AP[cong[0]].items():
        if len(m)==1 and m[0] in plift5.HANDLE: ck=abs(c)
    rows3.append((ca,cb,ck))
print("three conditions (c_k1, c_k2, modulus c_k):",rows3)

def resid(mu,ca,cb):
    m1,m2,m3,m4,m5,m6=mu
    a=m1+m2+m5; bb=m1-m2; d=m4-m3; g=m3+m6; h2=m2-m5
    dn1=(2*E*A*bb + a*A*A - 2*Bv*d) + P*(E*bb*bb + 2*a*A*bb - d*d) + P*P*a*bb*bb
    dn2=(A*g + H*bb - A*h2 - J*d) + P*(bb*g - d*h2)
    return ca*(n1+dn1)+cb*(n2+dn2)

print()
allq=set()
for ca,cb,ck in rows3: allq|=set(fac(ck))
print("primes over all three moduli:",sorted(allq))
overall=True
for ca,cb,ck in rows3:
    if ck==1:
        print("  modulus 1 -> condition vacuous"); continue
    for q,e in sorted(fac(ck).items()):
        qq=q**e
        found=None
        # single-parameter search
        for p_ in range(6):
            for t in range(qq):
                mu=[0]*6; mu[p_]=t*step[p_]
                if resid(mu,ca,cb)%qq==0: found=('1-param',p_,t); break
            if found: break
        if not found:
            for p1,p2 in combinations(range(6),2):
                for t1 in range(qq):
                    for t2 in range(qq):
                        mu=[0]*6; mu[p1]=t1*step[p1]; mu[p2]=t2*step[p2]
                        if resid(mu,ca,cb)%qq==0: found=('2-param',(p1,p2),(t1,t2)); break
                    if found: break
                if found: break
        print("  modulus %d, prime power %d -> solvable: %-5s %s"%(ck,qq,found is not None,found))
        if not found: overall=False
print()
print("BLOCK %d, |S|=2, 6-PARAMETER MODEL: all conditions simultaneously satisfiable per-prime:"%jj, overall)
