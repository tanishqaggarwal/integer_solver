import os, sys, itertools, collections
from fractions import Fraction
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad
P = 2**256-2**32-977
SEVEN = [22229, 22230, 35758, 35759, 35760, 35761, 35762]
E = [2554, 6816, 8124, 9123, 9421, 12231, 12270, 12350, 14584, 18673, 22044, 29125]
v = L.load(os.path.join(LAB,'best','new_instance_partial_39026.json'))
av = L.all_atom_values(v)
A = [av[a] for a in SEVEN]
C0 = (A[0] + 7376877*A[6]) % P
A10 = A[1] % P
print('C0 =', C0)
print('A1_0 =', A10, ' zero?', A10==0)

def rowsfor(atoms):
    R=[]
    for e in E:
        m,sq,co = L.eq_atoms[e]
        R.append([co.get(a,0) for a in atoms])
    return R

def kernel_Q(M, n):
    """rational kernel basis (integer primitive rows) of |M| x n matrix"""
    A_=[[Fraction(x) for x in r] for r in M]; nn=len(A_); piv=[]; r_=0
    for j in range(n):
        k=next((i for i in range(r_,nn) if A_[i][j]!=0), None)
        if k is None: continue
        A_[r_],A_[k]=A_[k],A_[r_]
        pv=A_[r_][j]; A_[r_]=[x/pv for x in A_[r_]]
        for i in range(nn):
            if i!=r_ and A_[i][j]!=0:
                f=A_[i][j]; A_[i]=[x-f*y for x,y in zip(A_[i],A_[r_])]
        piv.append(j); r_+=1
    ker=[]; ps=set(piv)
    import math
    for fc in range(n):
        if fc in ps: continue
        w=[Fraction(0)]*n; w[fc]=Fraction(1)
        for i,pj in enumerate(piv): w[pj]=-A_[i][fc]
        den=1
        for x in w: den = den*x.denominator//math.gcd(den,x.denominator)
        ker.append([int(x*den) for x in w])
    return ker

# ---------- n = 7 : all 6-subsets ----------
R7 = rowsfor(SEVEN)
print('\n=== n=7, all 6-subsets: does the 1-dim kernel meet BOTH congruences? ===')
hits=[]
cnt=collections.Counter()
for sel in itertools.combinations(range(12),6):
    K = kernel_Q([R7[i] for i in sel], 7)
    cnt[len(K)]+=1
    if len(K)!=1: 
        if len(K)>1: hits.append(('DIM>1',sel,K))
        continue
    w=K[0]
    u = (w[0] + 7376877*w[6]) % P    # coefficient of lambda in congruence 1
    t = w[1] % P                     # coefficient of lambda in congruence 2
    # need lambda: u*lam == C0, t*lam == A10  (mod p)
    ok=None
    if t==0 and u==0:
        ok = (C0%P==0 and A10%P==0)
    elif t==0:
        ok = (A10%P==0) and True   # lam = C0/u
    elif u==0:
        ok = (C0%P==0)
    else:
        ok = (C0*t - A10*u) % P == 0
    if ok: hits.append(('COMPAT',sel,w))
print('kernel-dim distribution over 6-subsets:', dict(cnt))
print('compatible 6-subsets:', len(hits))
for h in hits[:10]: print('   ', h[0], [E[i] for i in h[1]], h[2])

# also: 5-subsets, confirm dim 2 and that both congruences are met (sanity: witness)
print('\n=== n=7, 5-subsets: dim distribution + congruence feasibility ===')
cnt5=collections.Counter(); feas5=0; feas_examples=[]
for sel in itertools.combinations(range(12),5):
    K = kernel_Q([R7[i] for i in sel],7); cnt5[len(K)]+=1
    if len(K)<2: continue
    # solve  sum lam_j*(w_j0+7376877*w_j6) == C0 ;  sum lam_j*w_j1 == A10  (mod p)
    rows=[[ (w[0]+7376877*w[6])%P for w in K],[ w[1]%P for w in K]]
    rhs=[C0%P, A10%P]
    # gaussian elim mod p on 2 x len(K)
    m=len(K); M=[rows[0]+[rhs[0]], rows[1]+[rhs[1]]]; r_=0; 
    for j in range(m):
        k=next((i for i in range(r_,2) if M[i][j]%P), None)
        if k is None: continue
        M[r_],M[k]=M[k],M[r_]
        inv=pow(M[r_][j],-1,P); M[r_]=[x*inv%P for x in M[r_]]
        for i in range(2):
            if i!=r_ and M[i][j]%P:
                f=M[i][j]; M[i]=[(a-f*b)%P for a,b in zip(M[i],M[r_])]
        r_+=1
    bad=any(all(M[i][j]%P==0 for j in range(m)) and M[i][m]%P for i in range(2))
    if not bad: feas5+=1; feas_examples.append(sel)
print('kernel-dim distribution over 5-subsets:', dict(cnt5))
print('congruence-feasible 5-subsets:', feas5, 'of', sum(cnt5.values()))
print('  witness subset {2554,6816,8124,9123,9421} feasible?', tuple(range(5)) in set(feas_examples))
print('  examples:', [[E[i] for i in s] for s in feas_examples[:6]])
