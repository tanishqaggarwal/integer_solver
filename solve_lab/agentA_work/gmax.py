"""GENERIC regional max-satisfy: low-support search (mod q) + exact HNF integrality."""
import sys, json, random, time, collections
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentA_work')
import env, lib as L
from regsolve2 import build
P=env.P; Q=(1<<61)-1

def int_solve(rowsM,rhs,nk):
    nr=len(rowsM)
    H=[r[:] for r in rowsM]; U=[[1 if i==j else 0 for j in range(nk)] for i in range(nk)]
    piv=[]; r=0
    for i in range(nr):
        if r>=nk: break
        while True:
            nzc=[j for j in range(r,nk) if H[i][j]]
            if len(nzc)<=1: break
            nzc.sort(key=lambda j: abs(H[i][j])); j0=nzc[0]
            for j in nzc[1:]:
                q=H[i][j]//H[i][j0]
                if q:
                    for k in range(nr): H[k][j]-=q*H[k][j0]
                    for k in range(nk): U[k][j]-=q*U[k][j0]
        nzc=[j for j in range(r,nk) if H[i][j]]
        if not nzc: continue
        j0=nzc[0]
        if j0!=r:
            for k in range(nr): H[k][r],H[k][j0]=H[k][j0],H[k][r]
            for k in range(nk): U[k][r],U[k][j0]=U[k][j0],U[k][r]
        piv.append((i,r)); r+=1
    y=[0]*nk
    for i,j in piv:
        s=rhs[i]-sum(H[i][k]*y[k] for k in range(j))
        if s%H[i][j]: return None
        y[j]=s//H[i][j]
    for i in range(nr):
        if sum(H[i][k]*y[k] for k in range(nk))!=rhs[i]: return None
    return [sum(U[k][j]*y[j] for j in range(nk)) for k in range(nk)]

def run(path, tlimit=300, target=None, seed=1):
    v=L.load(path); av=L.all_atom_values(v); fe=L.failing_eqs(av)
    s0=L.NEQ-len(fe)
    A=set(a for e in fe for a in L.eq_atoms[e][2])
    K,R,rows=build(v,A); nk=len(K)
    good=[(e,c,lin) for e,c,lin,hq in rows if not hq]
    dead=[(e,c) for e,c,lin,hq in rows if not hq and not lin and c!=0]
    NZ=[(e,c,lin) for e,c,lin in good if lin]
    n=len(NZ); EQ=[e for e,_,_ in NZ]
    N=[[lin.get(j,0) for j in range(nk)] for e,c,lin in NZ]
    B=[-c for e,c,lin in NZ]
    Nq=[[x%Q for x in r] for r in N]
    d0=[v[u] for u in K]
    cur=[EQ[i] for i in range(n) if sum(N[i][j]*d0[j] for j in range(nk))!=B[i]]
    print('%s score=%d failing=%d | region knobs=%d eqs=%d nontrivial=%d dead-violated=%d in-region-violated=%d'%(
        path.split('/')[-1],s0,len(fe),nk,len(R),n,len(dead),len(cur)),flush=True)
    if not n or nk==0 or n<nk: 
        print('   (degenerate region, skip)'); return
    def infoset_codewords(I):
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
    best={}
    random.seed(seed); t0=time.time(); trials=0
    IDX=list(range(n)); lim = target if target else max(1,len(cur)-1)
    while time.time()-t0<tlimit:
        trials+=1
        I=random.sample(IDX,nk)
        cws=infoset_codewords(I)
        if cws is None: continue
        OUT=[i for i in IDX if i not in I]
        for t in range(nk):
            c=cws[t]; sup=frozenset(EQ[i] for i in IDX if c[i])
            if len(sup)<=lim+2: best.setdefault(sup,1)
        for t1 in range(nk):
            for t2 in range(t1+1,nk):
                c1,c2=cws[t1],cws[t2]
                for k in OUT:
                    if not c2[k]: continue
                    lam=(-c1[k])*pow(c2[k],Q-2,Q)%Q
                    sup=frozenset(EQ[i] for i in IDX if (c1[i]+lam*c2[i])%Q)
                    if 0<len(sup)<=lim+2: best.setdefault(sup,2)
    print('   ISD trials=%d supports=%d sizes=%s'%(trials,len(best),
          sorted(collections.Counter(len(s) for s in best).items())),flush=True)
    bestres=(len(cur),None,None)
    for sup in sorted(best,key=len):
        if len(sup)>=bestres[0]: break
        D=set(sup); Z=[i for i in range(n) if EQ[i] not in D]
        x=int_solve([N[i] for i in Z],[B[i] for i in Z],nk)
        print('     |D|=%d %s -> %s'%(len(sup),sorted(sup),'INTEGRAL' if x else 'no'),flush=True)
        if x is not None: bestres=(len(sup),sorted(sup),x)
    if bestres[2] is not None:
        w=list(v)
        for j,u in enumerate(K): w[u]=bestres[2][j]
        av2=L.all_atom_values(w); fe2=L.failing_eqs(av2); s2=L.NEQ-len(fe2)
        print('   *** region violated %d -> %d ; SCORE %d -> %d'%(len(cur),bestres[0],s0,s2),flush=True)
        if s2>s0:
            out='/home/user/integer_solver/solve_lab/agentA_work/A_gmax_%d.json'%s2
            json.dump({str(i):str(w[i]) for i in range(L.NVARS)},open(out,'w')); print('   saved',out,flush=True)
    else:
        print('   NO improvement found; %d violated stands'%len(cur),flush=True)

if __name__=='__main__':
    tl=float(sys.argv[1])
    for p in sys.argv[2:]: run(p,tlimit=tl)
