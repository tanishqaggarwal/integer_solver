"""For each boolean frame: exact reduce + solve the residual system over F_p.
Reports SOLVABLE / verdict."""
import os, sys, itertools, pickle
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentG_work')
import g29_frame as F
import gsym2 as G
from gsym2 import L, ad, P

def rref(M,nc):
    M=[r[:] for r in M]; piv=[]; r=0
    for c in range(nc):
        pr=None
        for i in range(r,len(M)):
            if M[i][c]%P: pr=i;break
        if pr is None: continue
        M[r],M[pr]=M[pr],M[r]
        iv=pow(M[r][c],-1,P); M[r]=[x*iv%P for x in M[r]]
        for i in range(len(M)):
            if i!=r and M[i][c]%P:
                f=M[i][c]; M[i]=[(x-f*y)%P for x,y in zip(M[i],M[r])]
        piv.append(c); r+=1
        if r==len(M): break
    return M,piv,r

def verdict(v):
    r=F.analyse(v)
    if r['incchecks']: return 'INCONSISTENT-LINEAR(%d) %s'%(len(r['incchecks']),[(a,len(L.atom2eq.get(a,{}))) for a in r['incchecks']][:5]), r
    if r['nzc']: return 'UNREACHABLE-CONST(%d) %s'%(len(r['nzc']),[(a,len(L.atom2eq.get(a,{}))) for a,_ in r['nzc']][:5]), r
    res=r['res']
    if not res: return '*** FULL MOD-P SOLUTION ***', r
    consts=[(a,g%P) for a,g in res if isinstance(g,int) and g%P]
    if consts: return 'CONST-RESIDUAL(%d)'%len(consts), r
    polys=[(a,g) for a,g in res if not isinstance(g,int)]
    vars_=sorted({F.NB[k] for a,g in polys for m in g for k,_ in m})
    ix={u:i for i,u in enumerate(vars_)}; nv=len(vars_)
    lin=[(a,g) for a,g in polys if G.deg(g)==1]
    M=[]
    for a,g in lin:
        row=[0]*(nv+1)
        for m,c in g.items():
            if not m: row[nv]=(-c)%P
            else: row[ix[F.NB[m[0][0]]]]=c%P
        M.append(row)
    if M:
        MM,piv,rk=rref(M,nv)
        if any(all(x%P==0 for x in MM[i][:nv]) and MM[i][nv]%P for i in range(len(MM))):
            return 'RESIDUAL-LIN-INCONSISTENT nv=%d'%nv, r
    else:
        piv=[]; rk=0
    if rk==nv:
        sol=[0]*nv
        for i,c in enumerate(piv): sol[c]=MM[i][nv]%P
        bad=[]
        for a,g in polys:
            val=0
            for m,c in g.items():
                t=c
                for k,e in m: t=t*pow(sol[ix[F.NB[k]]],e,P)%P
                val=(val+t)%P
            if val: bad.append(a)
        if not bad: return '*** SOLVABLE (unique point) ***', r
        return 'PINNED-THEN-FAILS(%d) vars=%s'%(len(bad),vars_), r
    return 'UNDERDETERMINED rank %d of %d vars %s : %d polys'%(rk,nv,vars_,len(polys)), r

if __name__=='__main__':
    CTRL=[int(x) for x in sys.argv[1].split(',')]
    MAXK=int(sys.argv[2]) if len(sys.argv)>2 else len(CTRL)
    for k in range(0,MAXK+1):
        for combo in itertools.combinations(CTRL,k):
            v=list(F.v0)
            for b in combo: v[b]=1-v[b]
            ad.fwd(v,rounds=8)
            try: vd,r=verdict(v)
            except Exception as e: vd='ERR %s'%str(e)[:60]
            print('%-32s %s'%(','.join(map(str,combo)) or '(base)', vd), flush=True)
