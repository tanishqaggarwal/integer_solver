"""Equation-level (off-manifold) windows: full structure per level.
Rows are EQUATION VALUES as exact affine forms in the knobs; atoms are free to be nonzero
and to cancel.  Report rank, Q-consistency, currently-violated count, and whether the
uniqueness lemma applies."""
import sys, collections, json, time; sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentA_work')
import env, lib as L
from regsolve2 import build
P=env.P; Q=(1<<61)-1
path=sys.argv[1]; LEV=int(sys.argv[2])
v=L.load(path); av=L.all_atom_values(v); fe=L.failing_eqs(av)
print('%s score=%d failing=%d'%(path.split('/')[-1],L.NEQ-len(fe),len(fe)),flush=True)
def rank_mod(rowsM,nk,aug=None,q=Q):
    M=[[x%q for x in r]+([aug[i]%q] if aug else []) for i,r in enumerate(rowsM)]
    nr=len(M); r=0
    for c in range(nk):
        pr=None
        for i in range(r,nr):
            if M[i][c]: pr=i;break
        if pr is None: continue
        M[r],M[pr]=M[pr],M[r]
        inv=pow(M[r][c],q-2,q); M[r]=[x*inv%q for x in M[r]]
        for i in range(nr):
            if i!=r and M[i][c]:
                f=M[i][c]; M[i]=[(a-f*b)%q for a,b in zip(M[i],M[r])]
        r+=1
    inc=sum(1 for i in range(r,nr) if aug and M[i][nk]) if aug else 0
    return r,inc
A=set(a for e in fe for a in L.eq_atoms[e][2])
out=[]
for lev in range(LEV+1):
    R=sorted(set(e for a in A for e in L.atom2eq[a]))
    K,Rr,rows=build(v,A); nk=len(K)
    good=[(e,c,lin) for e,c,lin,hq in rows if not hq]
    skipped=len(rows)-len(good)
    NZ=[(e,c,lin) for e,c,lin in good if lin]
    n=len(NZ); EQ=[e for e,_,_ in NZ]
    N=[[lin.get(j,0) for j in range(nk)] for e,c,lin in NZ]
    B=[-c for e,c,lin in NZ]
    d0=[v[u] for u in K]
    cur=[EQ[i] for i in range(n) if sum(N[i][j]*d0[j] for j in range(nk))!=B[i]]
    r,inc=rank_mod(N,nk,B)
    print('lev%d atoms=%-5d eqs=%-5d skipped=%-3d nontrivial=%-5d knobs=%-4d rank=%-4d Qincons=%-3d violated=%-3d lemma_applies=%s'%(
        lev,len(A),len(Rr),skipped,n,nk,r,inc,len(cur),(r==nk and inc==0)),flush=True)
    out.append({'lev':lev,'atoms':len(A),'eqs':len(Rr),'n':n,'nk':nk,'rank':r,'inc':inc,'viol':len(cur)})
    if lev<LEV:
        A=set(a for e in R for a in L.eq_atoms[e][2])
json.dump(out,open('/home/user/integer_solver/solve_lab/agentA_work/eqwin.json','w'))
