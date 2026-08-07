"""Region model v2: knobs restricted so that EVERY region atom is affine in the knobs
(at most one knob per monomial).  Guarantees a fully linear, exact model."""
import sys, json, collections; sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentA_work')
from fractions import Fraction as F
import env, lib as L
P=env.P

def pick_knobs(v, A, prefer_free=True):
    A=set(A)
    cand=set(u for u in set(u for a in A for u in L.avars[a]) if all(x in A for x in L.var_atoms[u]))
    # drop knobs until every monomial of every region atom has <=1 knob
    while True:
        conflict=collections.Counter(); found=False
        for a in A:
            for m in L.polys[a]:
                ks=[u for u in set(m) if u in cand]
                if len(ks)>1 or (len(ks)==1 and m.count(ks[0])>1):
                    found=True
                    for u in ks: conflict[u]+=1
        if not found: break
        # drop the most conflicting; tie-break: prefer keeping FREE inputs and small-atom vars
        best=max(conflict.items(), key=lambda kv:(kv[1], len(L.var_atoms[kv[0]]), 1 if v[kv[0]]==P else 0))
        cand.discard(best[0])
    return sorted(cand)

def build(v, A, K=None):
    A=sorted(set(A)); Aset=set(A)
    av=L.all_atom_values(v)
    if K is None: K=pick_knobs(v,A)
    ki={u:i for i,u in enumerate(K)}
    def aff(a):
        c=0; lin=collections.defaultdict(int)
        for m,cc in L.polys[a].items():
            ks=[u for u in m if u in ki]
            assert len(ks)<=1, (a,m)
            if not ks:
                t=cc
                for u in m: t*=v[u]
                c+=t
            else:
                t=cc
                for u in m:
                    if u!=ks[0]: t*=v[u]
                lin[ki[ks[0]]]+=t
        return c,dict(lin)
    AFF={a:aff(a) for a in A}
    R=sorted(set(e for a in A for e in L.atom2eq[a]))
    rows=[]
    for e in R:
        mm,sq,co=L.eq_atoms[e]
        c=0; lin=collections.defaultdict(int); ok=True
        for a,cc in co.items():
            if a not in Aset:
                if av[a]!=0: ok=False
                continue
            c0,l0=AFF[a]; c+=cc*c0
            for i,x in l0.items(): lin[i]+=cc*x
        rows.append((e,c,dict(lin),not ok))
    return K,R,rows

def qsolve(rows,nk):
    a=[(e,c,lin) for e,c,lin,hq in rows if not hq]
    mat=[[F(lin.get(j,0)) for j in range(nk)]+[F(-c)] for e,c,lin in a]
    nr=len(mat); r=0; piv=[]
    for col in range(nk):
        pr=None
        for i in range(r,nr):
            if mat[i][col]!=0: pr=i;break
        if pr is None: continue
        mat[r],mat[pr]=mat[pr],mat[r]
        pv=mat[r][col]; mat[r]=[x/pv for x in mat[r]]
        for i in range(nr):
            if i!=r and mat[i][col]!=0:
                f=mat[i][col]; mat[i]=[x-f*y for x,y in zip(mat[i],mat[r])]
        piv.append(col); r+=1
    incons=[a[i][0] for i in range(r,nr) if mat[i][nk]!=0]
    sol=[None]*nk
    for i,c in enumerate(piv): sol[c]=mat[i][nk]
    return sol,[j for j in range(nk) if sol[j] is None],incons,r,piv,mat,[x[0] for x in a]
