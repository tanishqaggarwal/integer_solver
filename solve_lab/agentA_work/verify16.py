"""Verify the Lemma's hypotheses at L=16 explicitly: rank(N) = #knobs, Q-consistency,
and that the unique rational solution is NOT integral."""
import sys, collections; sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentA_work')
from fractions import Fraction as F
import env, lib as L
from regsolve2 import build
P=env.P; Q=(1<<61)-1
v=L.load('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json')
fe=L.failing_eqs(L.all_atom_values(v))
A=set(a for e in fe for a in L.eq_atoms[e][2])
LEV=int(sys.argv[1])
for _ in range(LEV):
    R=set()
    for a in A: R|=set(L.atom2eq[a])
    A=set(a for e in R for a in L.eq_atoms[e][2])
K,Rr,rows=build(v,A); nk=len(K)
good=[(e,c,lin) for e,c,lin,hq in rows if not hq]
assert len(good)==len(rows), 'window not exactly affine'
NZ=[(e,c,lin) for e,c,lin in good if lin]
n=len(NZ)
N=[[lin.get(j,0)%Q for j in range(nk)] for e,c,lin in NZ]
B=[(-c)%Q for e,c,lin in NZ]
M=[N[i][:]+[B[i]] for i in range(n)]
r=0; piv=[]
for c in range(nk):
    pr=None
    for i in range(r,n):
        if M[i][c]: pr=i;break
    if pr is None: continue
    M[r],M[pr]=M[pr],M[r]
    inv=pow(M[r][c],Q-2,Q); M[r]=[x*inv%Q for x in M[r]]
    for i in range(n):
        if i!=r and M[i][c]:
            f=M[i][c]; M[i]=[(a-f*b)%Q for a,b in zip(M[i],M[r])]
    piv.append(c); r+=1
inc=sum(1 for i in range(r,n) if M[i][nk])
print('L=%d: atoms=%d eqs=%d nontrivial rows=%d knobs=%d'%(LEV,len(A),len(Rr),n,nk))
print('   rank(N) = %d  (== #knobs ? %s)'%(r,r==nk))
print('   rows inconsistent over Q: %d  (Q-consistent ? %s)'%(inc,inc==0))
# unique solution: is it integral?  solve exactly with Fractions on the pivot rows only
Nz=[[F(lin.get(j,0)) for j in range(nk)] for e,c,lin in NZ]
Bz=[F(-c) for e,c,lin in NZ]
aug=[Nz[i]+[Bz[i]] for i in range(n)]
r2=0; piv2=[]
for c in range(nk):
    pr=None
    for i in range(r2,n):
        if aug[i][c]!=0: pr=i;break
    if pr is None: continue
    aug[r2],aug[pr]=aug[pr],aug[r2]
    pv=aug[r2][c]; aug[r2]=[x/pv for x in aug[r2]]
    for i in range(n):
        if i!=r2 and aug[i][c]!=0:
            f=aug[i][c]; aug[i]=[x-f*y for x,y in zip(aug[i],aug[r2])]
    piv2.append(c); r2+=1
sol=[None]*nk
for i,c in enumerate(piv2): sol[c]=aug[i][nk]
bad=[(K[j],sol[j].denominator) for j in range(nk) if sol[j] is not None and sol[j].denominator!=1]
print('   exact rank over Q = %d ; free coords = %d'%(r2,sum(1 for s in sol if s is None)))
print('   unique solution W is integral ? %s ; non-integral coords = %d'%(not bad,len(bad)))
for u,d in bad[:6]:
    print('        x%-7d denominator = %s'%(u,'p' if d==P else ('p*%d'%(d//P) if d%P==0 else d)))
