"""EXACT test: can 6 of the 12 equations be satisfied simultaneously?

Achievable atom set (derived and verified below):
  A = (A0..A6) for atoms 22229,22230,35758,35759,35760,35761,35762
  A0 + 7376877*A6 == Ca  (mod 7376877*P)      Ca = x_7068 - x_2099
  A1            == Cb    (mod P)              Cb = x_28730
  A2..A5 free integers
For each 6-subset T of the 12 equations, the integer kernel of M_T is Z*w (w primitive).
A = t*w must satisfy the two congruences -> a pair of congruences in one unknown t.
"""
import os, sys, json, itertools
from fractions import Fraction
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L
P=2**256-2**32-977
S=[22229,22230,35758,35759,35760,35761,35762]
E0=sorted(set().union(*[set(L.atom2eq[a]) for a in S]))
v=L.load(os.path.join(LAB,'best','new_instance_partial_39026.json'))
av=L.all_atom_values(v)
A=[av[a] for a in S]
Ca=v[7068]-v[2099]
Cb=v[28730]
M=7376877*P
print('check invariants on delivered witness:')
print('  (A0+7376877*A6 - Ca) %% (7376877*P) =',(A[0]+7376877*A[6]-Ca)%M)
print('  (A1 - Cb) %% P =',(A[1]-Cb)%P)
print('  Ca %% P != 0 :',Ca%P!=0,'  Cb %% P != 0 :',Cb%P!=0)
sat=[e for e in E0 if L.eq_value(e,av)==0]
print('  currently satisfied:',sat,'(%d of %d)'%(len(sat),len(E0)))

def coefrow(e):
    m,sq,co=L.eq_atoms[e]
    return [co.get(a,0) for a in S]
ROWS={e:coefrow(e) for e in E0}

def kernel(rows):
    """integer primitive basis of the rational kernel of the matrix `rows` (list of lists)."""
    n=7
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
    freec=[c for c in range(n) if c not in piv]
    basis=[]
    for fc in freec:
        vec=[Fraction(0)]*n; vec[fc]=Fraction(1)
        for c,rr in piv.items(): vec[c]=-Mx[rr][fc]
        den=1
        for x in vec: den=den*x.denominator//__import__('math').gcd(den,x.denominator)
        iv=[int(x*den) for x in vec]
        g=0
        for x in iv: g=__import__('math').gcd(g,x)
        if g: iv=[x//g for x in iv]
        basis.append(iv)
    return basis

def solve_two(a1,m1,b1, a2,m2,b2):
    """t with a1*t = b1 (mod m1) and a2*t = b2 (mod m2); return t or None."""
    import math
    def lin(a,b,m):
        g=math.gcd(a%m,m)
        if b%g: return None
        a2_,b2_,m2_=(a%m)//g,(b%m)//g,m//g
        return (b2_*pow(a2_,-1,m2_))%m2_, m2_
    r1=lin(a1,b1,m1)
    if r1 is None: return None
    r2=lin(a2,b2,m2)
    if r2 is None: return None
    t1,n1=r1; t2,n2=r2
    g=math.gcd(n1,n2)
    if (t2-t1)%g: return None
    lcm=n1//g*n2
    t=(t1+ n1*(((t2-t1)//g*pow(n1//g,-1,n2//g))%(n2//g)))%lcm if n2//g>1 else t1%lcm
    return t%lcm

hits=[]
for T in itertools.combinations(E0,6):
    B=kernel([ROWS[e] for e in T])
    if len(B)!=1:
        hits.append(('DIM',T,len(B))); continue
    w=B[0]
    a1=w[0]+7376877*w[6]; a2=w[1]
    t=solve_two(a1,M,Ca, a2,P,Cb)
    if t is not None:
        hits.append(('HIT',T,w,t))
print('\n6-subsets: kernel-dim anomalies / hits:')
n_hit=0
for h in hits:
    if h[0]=='HIT':
        n_hit+=1
        print('  HIT',h[1],'w=',h[2],'t=',str(h[3])[:40])
    else:
        print('  kernel dim !=1:',h[1],h[2])
print('total compatible 6-subsets:',n_hit,'of 924')

# also 7-subsets for completeness (kernel dim 0 -> only A=0, impossible)
bad=0
for T in itertools.combinations(E0,7):
    B=kernel([ROWS[e] for e in T])
    if len(B)>0: bad+=1
print('7-subsets with nonzero kernel:',bad,'of',len(list(itertools.combinations(E0,7))))
json.dump([[h[0],list(h[1])] for h in hits],open(os.path.join(HERE,'pa_six.json'),'w'))
