"""EXACT test of the 8-atom placement (base 7 + a22231, paying a37887 = T^2).

Achievable set, derived from the circuit and verified numerically below:
  A = (A22229, A22230, A35758, A35759, A35760, A35761, A35762, A22231) in Z^8
  (i)  A22229 + 7376877*A35762 == Ca  (mod 7376877*P)     Ca = x_7068 - x_2099
  (ii) A22230 + A22231         == C1  (mod P)             C1 = x_4432 - x_19964
  A35758..A35761 free.
If some 7-subset of the 12 equations admits an A in this set, the score becomes
12-7 = 5 failing among the twelve plus 1 for a37887  =>  6 total  =>  39,027.
"""
import os, sys, json, itertools, math
from fractions import Fraction
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L
P=2**256-2**32-977; MM=7376877*P
S=[22229,22230,35758,35759,35760,35761,35762,22231]
E0=sorted(set().union(*[set(L.atom2eq[a]) for a in S]))
v=L.load(os.path.join(LAB,'best','new_instance_partial_39026.json'))
Ca=v[7068]-v[2099]; C1=v[4432]-v[19964]
print('|E0| =',len(E0),E0)
print('Ca mod P =',Ca%P,'   C1 mod P =',C1%P)

def rows_of(T):
    out=[]
    for e in T:
        m,sq,co=L.eq_atoms[e]
        out.append([co.get(a,0) for a in S])
    return out

def kernel(rows,n):
    Mx=[[Fraction(x) for x in r] for r in rows]
    piv={}; r=0
    for c in range(n):
        k=None
        for i in range(r,len(Mx)):
            if Mx[i][c]: k=i; break
        if k is None: continue
        Mx[r],Mx[k]=Mx[k],Mx[r]
        f=Mx[r][c]; Mx[r]=[x/f for x in Mx[r]]
        for i in range(len(Mx)):
            if i!=r and Mx[i][c]:
                g=Mx[i][c]; Mx[i]=[a-g*b for a,b in zip(Mx[i],Mx[r])]
        piv[c]=r; r+=1
    basis=[]
    for fc in [c for c in range(n) if c not in piv]:
        vec=[Fraction(0)]*n; vec[fc]=Fraction(1)
        for c,rr in piv.items(): vec[c]=-Mx[rr][fc]
        den=1
        for x in vec: den=den*x.denominator//math.gcd(den,x.denominator)
        iv=[int(x*den) for x in vec]
        g=0
        for x in iv: g=math.gcd(g,x)
        if g: iv=[x//g for x in iv]
        basis.append(iv)
    return basis

def lin(a,b,m):
    g=math.gcd(a%m,m)
    if b%g: return None
    return ((b//g)*pow((a%m)//g,-1,m//g))%(m//g), m//g

def crt(t1,n1,t2,n2):
    g=math.gcd(n1,n2)
    if (t2-t1)%g: return None
    l=n1//g*n2
    return (t1+n1*(((t2-t1)//g)*pow(n1//g,-1,n2//g)%(n2//g)))%l if n2//g>1 else t1%l

dims=collections_dims={}
hits=[]
for k in (8,7):
    cnt={}
    for T in itertools.combinations(E0,k):
        B=kernel(rows_of(T),8)
        cnt[len(B)]=cnt.get(len(B),0)+1
        if k==7 and len(B)==1:
            w=B[0]
            a1=w[0]+7376877*w[6]; a2=w[1]+w[7]
            r1=lin(a1,Ca,MM); r2=lin(a2,C1,P)
            if r1 and r2:
                t=crt(r1[0],r1[1],r2[0],r2[1])
                if t is not None: hits.append((T,w,t))
        elif k==7 and len(B)>1:
            hits.append(('DIM',T,len(B)))
    print(f'{k}-subsets kernel-dim histogram: {sorted(cnt.items())}  (n={sum(cnt.values())})')
print('\ncompatible 7-subsets (=> 6 failing total, 39027):',len(hits))
for h in hits[:10]: print('  ',h)
json.dump([[list(h[0]) if not isinstance(h[0],str) else h[0]] for h in hits],open(os.path.join(HERE,'pa_eight.json'),'w'))
