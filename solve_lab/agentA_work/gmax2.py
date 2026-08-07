"""gmax with equation-closure region expansion."""
import sys, json, random, time, collections
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentA_work')
import env, lib as L
from regsolve2 import build
from gmax import int_solve
P=env.P; Q=(1<<61)-1

def region(v, fe, ecl):
    A=set(a for e in fe for a in L.eq_atoms[e][2])
    for _ in range(ecl):
        R=set(e for a in A for e in L.atom2eq[a])
        A2=set(a for e in R for a in L.eq_atoms[e][2])
        if A2==A: break
        A=A2
    return A

def run(path, ecl=1, tlimit=180, seed=5, lim=None):
    v=L.load(path); av=L.all_atom_values(v); fe=L.failing_eqs(av); s0=L.NEQ-len(fe)
    A=region(v,fe,ecl)
    K,R,rows=build(v,A); nk=len(K)
    good=[(e,c,lin) for e,c,lin,hq in rows if not hq]
    skipped=len(rows)-len(good)
    NZ=[(e,c,lin) for e,c,lin in good if lin]
    n=len(NZ); EQ=[e for e,_,_ in NZ]
    N=[[lin.get(j,0) for j in range(nk)] for e,c,lin in NZ]
    B=[-c for e,c,lin in NZ]
    Nq=[[x%Q for x in r] for r in N]
    d0=[v[u] for u in K]
    cur=[EQ[i] for i in range(n) if sum(N[i][j]*d0[j] for j in range(nk))!=B[i]]
    deadviol=[e for e,c,lin,hq in rows if not hq and not lin and c!=0]
    print('%s ecl=%d score=%d fail=%d | atoms=%d knobs=%d eqs=%d skipped=%d rows=%d deadviol=%d inregion-viol=%d'%(
        path.split('/')[-1],ecl,s0,len(fe),len(A),nk,len(R),skipped,n,len(deadviol),len(cur)),flush=True)
    if nk==0 or n<nk: print('   degenerate'); return
    if lim is None: lim=max(1,min(len(cur),len(fe))-1)
    def infoset(I):
        M=[[Nq[i][j] for j in range(nk)]+[1 if k==t else 0 for t in range(len(I))] for k,i in enumerate(I)]
        r=0
        for c in range(nk):
            pr=None
            for i in range(r,len(M)):
                if M[i][c]: pr=i;break
            if pr is None: continue
            M[r],M[pr]=M[pr],M[r]
            inv=pow(M[r][c],Q-2,Q); M[r]=[x*inv%Q for x in M[r]]
            for i in range(len(M)):
                if i!=r and M[i][c]:
                    f=M[i][c]; M[i]=[(a-f*b)%Q for a,b in zip(M[i],M[r])]
            r+=1
        if r!=nk: return None
        Uc=[[M[k][nk+t] for t in range(len(I))] for k in range(nk)]
        return [[sum(Nq[i][j]*Uc[j][t] for j in range(nk))%Q for i in range(n)] for t in range(nk)]
    best={}; random.seed(seed); t0=time.time(); trials=0; IDX=list(range(n))
    while time.time()-t0<tlimit:
        trials+=1
        I=random.sample(IDX,nk)
        cws=infoset(I)
        if cws is None: continue
        OUT=[i for i in IDX if i not in I]
        for t in range(nk):
            sup=frozenset(EQ[i] for i in IDX if cws[t][i])
            if 0<len(sup)<=lim: best.setdefault(sup,1)
        for t1 in range(nk):
            for t2 in range(t1+1,nk):
                c1,c2=cws[t1],cws[t2]
                for k in OUT:
                    if not c2[k]: continue
                    lam=(-c1[k])*pow(c2[k],Q-2,Q)%Q
                    sup=frozenset(EQ[i] for i in IDX if (c1[i]+lam*c2[i])%Q)
                    if 0<len(sup)<=lim: best.setdefault(sup,2)
    print('   trials=%d supports<=%d=%d sizes=%s'%(trials,lim,len(best),
          sorted(collections.Counter(len(s) for s in best).items())[:6]),flush=True)
    for sup in sorted(best,key=len):
        D=set(sup); Z=[i for i in range(n) if EQ[i] not in D]
        x=int_solve([N[i] for i in Z],[B[i] for i in Z],nk)
        if x is None: continue
        w=list(v)
        for j,u in enumerate(K): w[u]=x[j]
        s2=L.NEQ-len(L.failing_eqs(L.all_atom_values(w)))
        print('   INTEGRAL |D|=%d -> SCORE %d'%(len(sup),s2),flush=True)
        if s2>=39026:
            out='/home/user/integer_solver/solve_lab/agentA_work/A_gmax2_%d.json'%s2
            json.dump({str(i):str(w[i]) for i in range(L.NVARS)},open(out,'w')); print('   saved',out,flush=True)
        if s2>s0: return
    print('   no improvement',flush=True)

if __name__=='__main__':
    ecl=int(sys.argv[1]); tl=float(sys.argv[2])
    for p in sys.argv[3:]: run(p,ecl=ecl,tlimit=tl)
