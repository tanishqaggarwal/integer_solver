"""Enlarge the movable set: for a given extra-atom list, build the affine model,
report rank/consistency over Q (mod a big prime for speed) and integrality."""
import sys, json, random, collections; sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentA_work')
import env, lib as L
from amk_model import build, knobpoly, v0, av0
P=env.P
Q=(1<<61)-1   # prime 2^61-1

def model(extra, verbose=True):
    A,K,R=build(extra); Aset=set(A)
    kp={a:knobpoly(a,K,v0) for a in A}
    QUAD=[a for a in A if any(len(m)>1 for m in kp[a])]
    def affine(a):
        Pp=kp[a]; c=Pp.get((),0); lin={}
        for m,cc in Pp.items():
            if len(m)==1: lin[m[0]]=cc
        return c,lin
    rows=[]
    for e in R:
        mm,sq,co=L.eq_atoms[e]
        c=0; lin=collections.defaultdict(int); hq=False
        for a,cc in co.items():
            if a not in Aset:
                if av0[a]!=0: hq=True
                continue
            if a in QUAD: hq=True; continue
            c0,l0=affine(a); c+=cc*c0
            for i,x in l0.items(): lin[i]+=cc*x
        rows.append((e,c,dict(lin),hq))
    return A,K,R,rows,QUAD

def rank_consistency(rows,nk,q=Q):
    """returns (rank, consistent, rows_used) mod q"""
    aff=[(e,c%q,{i:x%q for i,x in lin.items()}) for e,c,lin,hq in rows if not hq]
    mat=[[lin.get(j,0) for j in range(nk)]+[(-c)%q] for e,c,lin in aff]
    nr=len(mat); r=0; piv=[]
    for col in range(nk):
        pr=None
        for i in range(r,nr):
            if mat[i][col]: pr=i;break
        if pr is None: continue
        mat[r],mat[pr]=mat[pr],mat[r]
        inv=pow(mat[r][col],q-2,q)
        mat[r]=[x*inv%q for x in mat[r]]
        for i in range(nr):
            if i!=r and mat[i][col]:
                f=mat[i][col]; mat[i]=[(a-f*b)%q for a,b in zip(mat[i],mat[r])]
        piv.append(col); r+=1
    incons=sum(1 for i in range(r,nr) if mat[i][nk])
    return r,len(aff),incons

if __name__=='__main__':
    import sys
    tests={'base':[], '+37887':[37887], '+37887,41906':[37887,41906],
           '+29090':[37887,41906,29090],
           '+29426':[37887,41906,29426],
           '+41972':[37887,41906,41972],
           '+36085':[37887,41906,36085],
           '+all4':[37887,41906,29426,41972],
           '+29090,29426,41972':[37887,41906,29090,29426,41972]}
    for name,ex in tests.items():
        A,K,R,rows,QUAD=model(ex)
        r,na,inc=rank_consistency(rows,len(K))
        nq=sum(1 for x in rows if x[3])
        print('%-24s atoms=%-4d knobs=%-4d eqs=%-5d affrows=%-5d quadrows=%-3d rank=%-4d Qincons=%d'%(
            name,len(A),len(K),len(R),na,nq,r,inc))
