"""Per-bit: solve the MAIN linear system AND the residual linear obligations, then read
the six coordinates, map to secp256k1, and test whether the freed coordinate lands on a
curve point (i.e. the bit selects a table point)."""
import os, sys, itertools, pickle, time
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentG_work')
import g29_frame as F, gpt
import gsym2 as G
from gsym2 import L, ad, P
NB=F.NB; n=len(NB); ixm={u:i for i,u in enumerate(NB)}

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

def frame(FL):
    v=list(F.v0)
    for b in FL: v[b]=1-v[b]
    ad.fwd(v,rounds=8)
    r=F.analyse(v)
    val,_=G.build(v,NB,cap=6)
    f1=G.evalatom(19297,val,6); f2=G.evalatom(19299,val,6)
    o={'flips':FL,'nres':len(r['res']),'ninc':len(r['incchecks']),'nzc':len(r['nzc'])}
    if isinstance(f1,int) or isinstance(f2,int): o['core']='DEAD'; return o
    Ap,Bp=gpt.pencil(f1,f2); lab=gpt.label(Ap,Bp,NB)
    if lab is None: o['core']='UNLABELLED'; return o
    o['label']=lab
    piv,R=r['piv'],r['R']
    t={c:(v[NB[c]]%P) for c in range(n) if c not in piv}
    # residual linear obligations: solve them for the still-free variables
    rlin=[(a,g) for a,g in r['res'] if not isinstance(g,int) and G.deg(g)==1]
    rvars=sorted({k for a,g in rlin for m in g for k,_ in m})
    if rvars:
        jx={c:i for i,c in enumerate(rvars)}; nv=len(rvars)
        M=[]
        for a,g in rlin:
            row=[0]*(nv+1)
            for m,c in g.items():
                if not m: row[nv]=(-c)%P
                else: row[jx[m[0][0]]]=c%P
            M.append(row)
        MM,rp,rk=rref(M,nv)
        o['res_lin_rank']=rk; o['res_lin_nvars']=nv
        o['res_lin_inconsistent']=any(all(x%P==0 for x in MM[i][:nv]) and MM[i][nv]%P for i in range(len(MM)))
        for i,c in enumerate(rp):
            if all(MM[i][cc]%P==0 for cc in range(nv) if cc!=c):
                t[rvars[c]]=MM[i][nv]%P
    def value(u):
        j=ixm.get(u)
        if j is None: return v[u]%P
        if j not in piv: return t.get(j,0)
        row=R[piv[j]]; val_=row.get(n,0)%P
        for c,vv in row.items():
            if c!=n and c!=j: val_=(val_-vv*t.get(c,0))%P
        return val_
    co={k:value(u) for k,u in lab.items()}
    o['coord']=co
    Q=[gpt.tosec(co['x1'],co['y1']),gpt.tosec(co['x2'],co['y2']),gpt.tosec(co['x3'],co['y3'])]
    o['pts']=Q
    o['oncurve']=[(q[1]*q[1]-pow(q[0],3,P)-7)%P==0 for q in Q]
    if all(o['oncurve']): o['D']=gpt.sub(Q[2],gpt.add(Q[0],Q[1]))
    return o

if __name__=='__main__':
    BITS=[int(x) for x in sys.argv[1].split(',')] if len(sys.argv)>1 else [47,91,112,438,490,542,853,1203,1357,1413,1438,1502]
    base=frame([]); print('BASE oncurve=%s'%base['oncurve']); print('  P1..P3:',base['pts']); print('  D=',base['D'],flush=True)
    res={():base}
    for b in BITS:
        try: o=frame([b])
        except Exception as e: print('x%d ERR %s'%(b,str(e)[:60])); continue
        res[(b,)]=o
        print('x%-6d nres=%-2d resLIN rank=%s/%s inc=%s oncurve=%s'%(b,o['nres'],o.get('res_lin_rank'),o.get('res_lin_nvars'),o.get('res_lin_inconsistent'),o.get('oncurve')),flush=True)
        if o.get('oncurve') and all(o['oncurve']):
            print('        P1=%s'%(o['pts'][0],))
            print('        P2=%s'%(o['pts'][1],))
            print('        D =%s'%(o.get('D'),),flush=True)
    pickle.dump(res,open('/home/user/integer_solver/solve_lab/agentG_work/tablescan.pkl','wb'))
