"""EXACT integer solvability of a state's residual region (HNF), not just the
denominators of one particular rational solution."""
import sys, json, collections, time; sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentA_work')
import env, lib as L
from regsolve2 import build, qsolve, pick_knobs
P=env.P

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
        if s%H[i][j]: return None,(i,j)
        y[j]=s//H[i][j]
    for i in range(nr):
        if sum(H[i][k]*y[k] for k in range(nk))!=rhs[i]: return None,('check',i)
    return [sum(U[k][j]*y[j] for j in range(nk)) for k in range(nk)],None

def run(path,grow=0):
    v=L.load(path); av=L.all_atom_values(v); fe=L.failing_eqs(av)
    print('%-38s score=%d failing=%d'%(path.split('/')[-1],L.NEQ-len(fe),len(fe)),flush=True)
    A=set(a for e in fe for a in L.eq_atoms[e][2])
    for g in range(grow+1):
        K,R,rows=build(v,A)
        nk=len(K)
        good=[(e,c,lin) for e,c,lin,hq in rows if not hq]
        skipped=len(rows)-len(good)
        M=[[lin.get(j,0) for j in range(nk)] for e,c,lin in good]
        B=[-c for e,c,lin in good]
        t0=time.time()
        x,why=int_solve(M,B,nk)
        print('   g%d atoms=%d knobs=%d eqs=%d skipped=%d -> %s  [%.1fs]'%(
            g,len(A),nk,len(R),skipped,'INTEGRAL SOLUTION' if x else 'NO integer solution (%s)'%str(why),time.time()-t0),flush=True)
        if x is not None:
            w=list(v)
            for j,u in enumerate(K): w[u]=x[j]
            av2=L.all_atom_values(w); fe2=L.failing_eqs(av2)
            print('   -> applied: score %d (was %d)'%(L.NEQ-len(fe2),L.NEQ-len(fe)),flush=True)
            s=L.NEQ-len(fe2)
            if s>=39026:
                out='/home/user/integer_solver/solve_lab/agentA_work/A_zsolve_%d.json'%s
                json.dump({str(i):str(w[i]) for i in range(L.NVARS)},open(out,'w'))
                print('   saved',out,flush=True)
            return w
        if g<grow:
            A |= set(x2 for a in list(A) for u in L.avars[a] for x2 in L.var_atoms[u])
    return None

if __name__=='__main__':
    for p in sys.argv[1:]: run(p,grow=0)
