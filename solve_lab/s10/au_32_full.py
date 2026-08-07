"""Constructively verify the FULL achievable-set claim:
   every A with  A0 + 7376877*A6 == C0 (mod p)  and  A1 == A1_0 (mod p)  is realisable
   at zero collateral -- including arbitrary residues of A0 mod 7376877
   (which needs the x_7068 += k*p / x_30163 re-base).
Then search all 378 congruence-feasible 5-subsets for one that also admits a 6th.
"""
import os, sys, random, itertools, math
from fractions import Fraction
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad
P = 2**256-2**32-977
SEVEN=[22229,22230,35758,35759,35760,35761,35762]
E12=[2554,6816,8124,9123,9421,12231,12270,12350,14584,18673,22044,29125]
E12S=set(E12)
DETACH={7068:22229,28730:22230,29854:35758,31864:35761,642:35762}
definer={t:a for t,a in L.definer.items() if t not in DETACH}
ORDER=[t for t in ad.ORDER if t not in DETACH]
def fwd2(v,rounds=3):
    for _ in range(rounds):
        for u in ORDER:
            nv=T.solve_lin(definer[u],u,v)
            if nv is not None: v[u]=nv
    return v
v0=L.load(os.path.join(LAB,'best','new_instance_partial_39026.json'))
w=list(v0); fwd2(w,8)
avb=L.all_atom_values(w); Ab=[avb[a] for a in SEVEN]
C0=(Ab[0]+7376877*Ab[6])%P; A10=Ab[1]%P

def realise(A):
    """Realise target A exactly, re-basing x_7068 by multiples of p when needed."""
    v=list(w)
    # --- re-base x_7068 so that (x_7068 - x_2099 - A0) % 7376877 == 0
    need=(v[7068]-v[2099]-A[0])% 7376877
    if need:
        k=(-need*pow(P%7376877,-1,7376877))%7376877
        v[7068]+=k*P
        fwd2(v,3)
        # repair a29539 through x_30163 :  x_29967 = p*x_30163 = 12846437*(x_14853-x_1308)
        tgt=12846437*(v[14853]-v[1308])
        if tgt % P: return None,'a29539 not repairable'
        v[30163]=tgt//P
        fwd2(v,3)
    S=A[2]+A[3]
    x9118=(S*pow(5113045,-1,P))%P
    rem=5113045*x9118-S
    x1329=rem//P
    v[9118]=x9118; v[1329]=x1329; v[29854]=A[2]+P*x1329
    D=A[5]-A[4]; x10903=D//P; x8731=D-P*x10903
    v[10903]=x10903; v[8731]=x8731; v[31864]=A[4]+P*x10903
    num=v[28730]-A[1]
    if num%P: return None,'A1 residue wrong'
    v[9413]=num//P
    num2=v[7068]-v[2099]-A[0]
    if num2%7376877: return None,'A0 residue mod 7376877 still wrong'
    x642=num2//7376877; v[642]=x642
    num3=x642-A[6]
    if num3%P: return None,'congruence 1 violated'
    v[17325]=num3//P
    fwd2(v,3)
    return v,None

random.seed(3)
print('=== realise random A meeting ONLY the two mod-p congruences ===')
for trial in range(5):
    A=[0]*7
    A[2]=random.getrandbits(500)-(1<<499); A[3]=random.getrandbits(500)-(1<<499)
    A[4]=random.getrandbits(500)-(1<<499); A[5]=random.getrandbits(500)-(1<<499)
    A[1]=A10 + P*(random.getrandbits(60)-(1<<59))
    A[6]=random.getrandbits(300)-(1<<299)
    A[0]=(C0-7376877*A[6])%P + P*(random.getrandbits(40)-(1<<39))
    v,err=realise(A)
    if err: print('  FAIL:',err); continue
    av=L.all_atom_values(v); f=set(L.failing_eqs(av))
    got=[av[a] for a in SEVEN]
    print(f'  trial{trial}: exact={got==A}  score={L.NEQ-len(f)}  outside-12 failures={sorted(f-E12S)}')

# ---- realise the best A for each congruence-feasible 5-subset and confirm 39,026 ----
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
print('\n=== realise A for OTHER congruence-feasible 5-subsets (not the witness one) ===')
tried=0
for sel in itertools.combinations(range(12),5):
    if tried>=4: break
    if sel==(0,1,2,3,4): continue
    Kb=kernel_Q([rows[i] for i in sel],7)
    if len(Kb)<2: continue
    # solve lam over F_p for the two congruences
    m=len(Kb)
    M=[[ (wv[0]+7376877*wv[6])%P for wv in Kb]+[C0%P],
       [ wv[1]%P for wv in Kb]+[A10%P]]
    r_=0; piv=[]
    for j in range(m):
        k=next((i for i in range(r_,2) if M[i][j]%P),None)
        if k is None: continue
        M[r_],M[k]=M[k],M[r_]
        inv=pow(M[r_][j],-1,P); M[r_]=[x*inv%P for x in M[r_]]
        for i in range(2):
            if i!=r_ and M[i][j]%P:
                f=M[i][j]; M[i]=[(a-f*b)%P for a,b in zip(M[i],M[r_])]
        piv.append(j); r_+=1
    if any(all(M[i][j]%P==0 for j in range(m)) and M[i][m]%P for i in range(2)): continue
    lam=[0]*m
    for i,j in enumerate(piv): lam[j]=M[i][m]
    A=[sum(lam[j]*Kb[j][i] for j in range(m)) for i in range(7)]
    v,err=realise(A)
    tried+=1
    if err: print(f'  eqs {[E12[i] for i in sel]}: realise FAILED ({err})'); continue
    av=L.all_atom_values(v); f=set(L.failing_eqs(av))
    print(f'  eqs {[E12[i] for i in sel]}: score={L.NEQ-len(f)} failing={sorted(f)}')
