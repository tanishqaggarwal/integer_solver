"""Robustness: rerun both exhaustive subset tests with the WEAKEST plausible congruence
model (mod P only, dropping the 7376877 factor).  Still 0 hits => the negative result
does not depend on how tightly the achievable set is modelled."""
import os, sys, itertools, math
from fractions import Fraction
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L
P=2**256-2**32-977
v=L.load(os.path.join(LAB,'best','new_instance_partial_39026.json'))
Ca=(v[7068]-v[2099])%P; C1=(v[4432]-v[19964])%P
def kern(rows,n):
    Mx=[[Fraction(x) for x in r] for r in rows]; piv={}; r=0
    for c in range(n):
        k=next((i for i in range(r,len(Mx)) if Mx[i][c]),None)
        if k is None: continue
        Mx[r],Mx[k]=Mx[k],Mx[r]; f=Mx[r][c]; Mx[r]=[x/f for x in Mx[r]]
        for i in range(len(Mx)):
            if i!=r and Mx[i][c]:
                g=Mx[i][c]; Mx[i]=[a-g*b for a,b in zip(Mx[i],Mx[r])]
        piv[c]=r; r+=1
    B=[]
    for fc in [c for c in range(n) if c not in piv]:
        vec=[Fraction(0)]*n; vec[fc]=Fraction(1)
        for c,rr in piv.items(): vec[c]=-Mx[rr][fc]
        d=1
        for x in vec: d=d*x.denominator//math.gcd(d,x.denominator)
        B.append([int(x*d) for x in vec])
    return B
for S,k,tag in ((( [22229,22230,35758,35759,35760,35761,35762]),6,'7 atoms / 6-subsets'),
                (([22229,22230,35758,35759,35760,35761,35762,22231]),7,'8 atoms / 7-subsets')):
    E0=sorted(set().union(*[set(L.atom2eq[a]) for a in S]))
    n=len(S); hits=0; tot=0; dims={}
    for T in itertools.combinations(E0,k):
        rows=[[L.eq_atoms[e][2].get(a,0) for a in S] for e in T]
        B=kern(rows,n); tot+=1; dims[len(B)]=dims.get(len(B),0)+1
        if len(B)!=1: hits+=1; continue          # dim>1 would be a hit
        w=B[0]
        a1=(w[0]+7376877*w[6])%P
        a2=(w[1]+(w[7] if n==8 else 0))%P
        if a1==0 or a2==0:
            if (a1==0 and Ca%P) or (a2==0 and C1%P): continue
            hits+=1; continue
        t1=Ca*pow(a1,-1,P)%P; t2=C1*pow(a2,-1,P)%P
        if t1==t2: hits+=1
    print(f'{tag}: kernel dims {sorted(dims.items())}, compatible-under-relaxed-model = {hits} of {tot}')
