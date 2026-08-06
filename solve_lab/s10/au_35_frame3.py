"""Frame 2+3: detach x_4432 as well.  Verify that moving x_28730 then breaks
a37887 ONLY (claim (4)'s cheapest escape), and construct the best state reachable:
a 6-subset of the twelve satisfied with B = a22231 != 0, giving 6 + 1 = 7 failures.
"""
import os, sys, itertools, math, random
from fractions import Fraction
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad
P = 2**256-2**32-977
SEVEN=[22229,22230,35758,35759,35760,35761,35762]; SS=set(SEVEN)
E12=[2554,6816,8124,9123,9421,12231,12270,12350,14584,18673,22044,29125]; E12S=set(E12)
DET3={7068:22229,28730:22230,29854:35758,31864:35761,642:35762,4432:22231}
definer={t:a for t,a in L.definer.items() if t not in DET3}
ORDER=[t for t in ad.ORDER if t not in DET3]
def fwd3(v,rounds=3):
    for _ in range(rounds):
        for u in ORDER:
            nv=T.solve_lin(definer[u],u,v)
            if nv is not None: v[u]=nv
    return v
v0=L.load(os.path.join(LAB,'best','new_instance_partial_39026.json'))
w=list(v0); fwd3(w,8)
av=L.all_atom_values(w); f=set(L.failing_eqs(av))
print('frame3 base: score',L.NEQ-len(f),'identical to witness:', w==v0)
for d in (1, 7, 12345, P, 3*P+5):
    w2=list(w); w2[28730]+=d; fwd3(w2)
    a2=L.all_atom_values(w2); f2=set(L.failing_eqs(a2))
    nzo=[a for a in range(L.NA) if a2[a] and a not in SS]
    print(f'  x_28730 += {d if d<10**6 else "p-ish"}: nonzero outside seven = {nzo}  '
          f'score={L.NEQ-len(f2)} out12={sorted(f2-E12S)}')
# --- best configuration: n=8, c=2, choose a 6-subset with kernel dim 2 meeting both congruences
A=[av[a] for a in SEVEN]
C0=(A[0]+7376877*A[6])%P
K=(w[4432]-w[19964])%P
rows=[]
for e in E12:
    m,sq,co=L.eq_atoms[e]; rows.append([co.get(a,0) for a in SEVEN+[22231]])
def kernel_Q(M,n):
    A_=[[Fraction(x) for x in r] for r in M]; nn=len(A_); piv=[]; r_=0
    for j in range(n):
        k=next((i for i in range(r_,nn) if A_[i][j]!=0),None)
        if k is None: continue
        A_[r_],A_[k]=A_[k],A_[r_]
        pv=A_[r_][j]; A_[r_]=[x/pv for x in A_[r_]]
        for i in range(nn):
            if i!=r_ and A_[i][j]!=0:
                fq=A_[i][j]; A_[i]=[x-fq*y for x,y in zip(A_[i],A_[r_])]
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
def realise8(tgt):
    """tgt = [A0..A6, B].  frame3: x_4432 fixed, x_28730 free -> B = K' - x_28730 exactly."""
    v=list(w)
    Kz = v[4432]-v[19964]
    v[28730]= Kz - tgt[7]
    num=v[28730]-tgt[1]
    if num%P: return None,'A1/B congruence violated'
    v[9413]=num//P
    S=tgt[2]+tgt[3]
    x9118=(S*pow(5113045,-1,P))%P; rem=5113045*x9118-S; x1329=rem//P
    v[9118]=x9118; v[1329]=x1329; v[29854]=tgt[2]+P*x1329
    D=tgt[5]-tgt[4]; x10903=D//P; x8731=D-P*x10903
    v[10903]=x10903; v[8731]=x8731; v[31864]=tgt[4]+P*x10903
    need=(v[7068]-v[2099]-tgt[0])%7376877
    if need:
        k=(-need*pow(P%7376877,-1,7376877))%7376877
        v[7068]+=k*P; fwd3(v,3)
        t2=12846437*(v[14853]-v[1308])
        if t2%P: return None,'a29539 unrepairable'
        v[30163]=t2//P; fwd3(v,3)
    num2=v[7068]-v[2099]-tgt[0]
    if num2%7376877: return None,'A0 mod 7376877'
    x642=num2//7376877; v[642]=x642
    num3=x642-tgt[6]
    if num3%P: return None,'congruence 1'
    v[17325]=num3//P
    fwd3(v,3)
    return v,None
best=None
for sel in itertools.combinations(range(12),6):
    Kb=kernel_Q([rows[i] for i in sel],8)
    if len(Kb)<2: continue
    M=[[ (wv[0]+7376877*wv[6])%P for wv in Kb]+[C0%P],
       [ (wv[1]+wv[7])%P for wv in Kb]+[K%P]]
    m=len(Kb); r_=0; piv=[]
    for j in range(m):
        k=next((i for i in range(r_,2) if M[i][j]%P),None)
        if k is None: continue
        M[r_],M[k]=M[k],M[r_]
        inv=pow(M[r_][j],-1,P); M[r_]=[x*inv%P for x in M[r_]]
        for i in range(2):
            if i!=r_ and M[i][j]%P:
                fq=M[i][j]; M[i]=[(a-fq*b)%P for a,b in zip(M[i],M[r_])]
        piv.append(j); r_+=1
    if any(all(M[i][j]%P==0 for j in range(m)) and M[i][m]%P for i in range(2)): continue
    lam=[0]*m
    for i,j in enumerate(piv): lam[j]=M[i][m]
    tgt=[sum(lam[j]*Kb[j][i] for j in range(m)) for i in range(8)]
    v,err=realise8(tgt)
    if err: continue
    a2=L.all_atom_values(v); f2=set(L.failing_eqs(a2))
    got=[a2[a] for a in SEVEN]+[a2[22231]]
    sc=L.NEQ-len(f2)
    if best is None or sc>best[0]:
        best=(sc,[E12[i] for i in sel],sorted(f2),got==tgt,list(v))
        print(f'  6-subset {[E12[i] for i in sel]}: score={sc} exact={got==tgt} failing={sorted(f2)}', flush=True)
    if sc>39026:
        T.save(v, os.path.join(HERE,'au_best.json'))
        print('*** SAVED au_best.json score',sc); break
print('\nBEST frame-3 result:', best[:4] if best else None)
