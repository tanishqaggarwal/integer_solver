"""Cheapest confirmation of the Lemma's third hypothesis at a level: the FULL row system
has no integer solution.  Sufficient: it is already inconsistent mod p (|D| = 0 fails)."""
import sys; sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentA_work')
import env, lib as L
from regsolve2 import build
P=env.P
v=L.load('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json')
fe=L.failing_eqs(L.all_atom_values(v))
for LEV in [int(x) for x in sys.argv[1:]]:
    A=set(a for e in fe for a in L.eq_atoms[e][2])
    for _ in range(LEV):
        R=set()
        for a in A: R|=set(L.atom2eq[a])
        A=set(a for e in R for a in L.eq_atoms[e][2])
    K,Rr,rows=build(v,A); nk=len(K)
    NZ=[(e,c,lin) for e,c,lin,hq in rows if not hq and lin]
    n=len(NZ)
    M=[[lin.get(j,0)%P for j in range(nk)]+[(-c)%P] for e,c,lin in NZ]
    r=0
    for c in range(nk):
        pr=None
        for i in range(r,n):
            if M[i][c]: pr=i;break
        if pr is None: continue
        M[r],M[pr]=M[pr],M[r]
        inv=pow(M[r][c],-1,P); M[r]=[x*inv%P for x in M[r]]
        for i in range(n):
            if i!=r and M[i][c]:
                f=M[i][c]; M[i]=[(a-f*b)%P for a,b in zip(M[i],M[r])]
        r+=1
    inc=sum(1 for i in range(r,n) if M[i][nk])
    print('L=%-3d rows=%-5d knobs=%-5d : full system inconsistent MOD P ? %s (%d bad rows)'%(
        LEV,n,nk,inc>0,inc))
    print('       => no integer solution of the full system => W is not integral. %s'%
          ('CONFIRMED' if inc>0 else 'NOT confirmed this way'))
