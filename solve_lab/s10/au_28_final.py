import os, sys, itertools, collections, math
from fractions import Fraction
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad
P = 2**256-2**32-977
SEVEN=[22229,22230,35758,35759,35760,35761,35762]
E12=[2554,6816,8124,9123,9421,12231,12270,12350,14584,18673,22044,29125]
E12S=frozenset(E12)
v=L.load(os.path.join(LAB,'best','new_instance_partial_39026.json'))
av=L.all_atom_values(v)
A=[av[a] for a in SEVEN]
C0=(A[0]+7376877*A[6])%P; A10=A[1]%P
rows=[]
for e in E12:
    m,sq,co=L.eq_atoms[e]; rows.append([co.get(a,0) for a in SEVEN])
def kernel_Q(M,n):
    A_=[[Fraction(x) for x in r] for r in M]; nn=len(A_); piv=[]; r_=0
    for j in range(n):
        k=next((i for i in range(r_,nn) if A_[i][j]!=0),None)
        if k is None: continue
        A_[r_],A_[k]=A_[k],A_[r_]
        pv=A_[r_][j]; A_[r_]=[x/pv for x in A_[r_]]
        for i in range(nn):
            if i!=r_ and A_[i][j]!=0:
                f=A_[i][j]; A_[i]=[x-f*y for x,y in zip(A_[i],A_[r_])]
        piv.append(j); r_+=1
    ker=[];ps=set(piv)
    for fc in range(n):
        if fc in ps: continue
        wv=[Fraction(0)]*n; wv[fc]=Fraction(1)
        for i,pj in enumerate(piv): wv[pj]=-A_[i][fc]
        den=1
        for x in wv: den=den*x.denominator//math.gcd(den,x.denominator)
        ker.append([int(x*den) for x in wv])
    return ker
def feasible(Kb, conds):
    m=len(Kb); M=[]
    for cvec,rhs in conds:
        M.append([sum(cvec[i]*wv[i] for i in range(len(wv)))%P for wv in Kb]+[rhs%P])
    nr=len(M); r_=0
    for j in range(m):
        k=next((i for i in range(r_,nr) if M[i][j]%P),None)
        if k is None: continue
        M[r_],M[k]=M[k],M[r_]
        inv=pow(M[r_][j],-1,P); M[r_]=[x*inv%P for x in M[r_]]
        for i in range(nr):
            if i!=r_ and M[i][j]%P:
                f=M[i][j]; M[i]=[(a-f*b)%P for a,b in zip(M[i],M[r_])]
        r_+=1
    for i in range(nr):
        if all(M[i][j]%P==0 for j in range(m)) and M[i][m]%P: return False
    return True
c1=[1,0,0,0,0,0,7376877]; c2=[0,1,0,0,0,0,0]
print('=== what does DROPPING a congruence actually buy? (n=7) ===')
for tag,conds in [('both (c=2)',[(c1,C0),(c2,A10)]),
                  ('only c1  (c=1)',[(c1,C0)]),
                  ('only c2  (c=1)',[(c2,A10)]),
                  ('none     (c=0)',[])]:
    best=None
    for k in range(12,0,-1):
        ok=0; ex=None
        for sel in itertools.combinations(range(12),k):
            Kb=kernel_Q([rows[i] for i in sel],7)
            if not Kb: continue
            if feasible(Kb,conds):
                ok+=1
                if ex is None: ex=sel
        if ok: best=(k,ok,ex); break
    k,ok,ex=best
    print(f'  {tag}: max satisfiable = {k} ({ok} subsets)  -> failing among the 12 = {12-k}   e.g. {[E12[i] for i in ex]}')
print()
print('=== compensators for the cheapest severing atoms ===')
for X in (21615, 30236, 32908, 32909, 38519, 32910, 21112, 3576):
    EX=frozenset(L.atom2eq[X])
    comp=[a for a in range(L.NA) if a!=X and frozenset(L.atom2eq[a])<=EX]
    print(f'  a{X}: neq={len(EX)} compensators-inside={comp}  ->  best-case |E|-n = {len(EX)-1-len(comp)}')
