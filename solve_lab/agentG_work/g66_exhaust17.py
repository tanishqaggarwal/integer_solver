"""EXHAUSTIVE minimum-cost departure on a support: group the affine rows into distinct
directions with multiplicities, then enumerate every violated-set of total multiplicity
<= BUDGET whose complement has smaller rank (so a departure exists), and test whether the
cubics can be zeroed on it."""
import os, sys, pickle, itertools, collections, time
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentG_work')
import gsym2 as G
from gsym2 import L, ad, P
import flint
D=pickle.load(open('/home/user/integer_solver/solve_lab/agentG_work/coset_model.pkl','rb'))
Lin=pickle.load(open('/home/user/integer_solver/solve_lab/agentG_work/coset_lin.pkl','rb'))
NB=D['NB']; lin=D['lin']; non=D['non']; x0=Lin['x0']; pt=D['pt']
ix={u:i for i,u in enumerate(NB)}
SUP=[ix[int(u)] for u in sys.argv[1].split(',')]
BUDGET=int(sys.argv[2]) if len(sys.argv)>2 else 7
k=len(SUP); sidx={c:j for j,c in enumerate(SUP)}
def sub(f):
    out={}
    for m,c in f.items():
        t=c%P; e=[0]*k
        for col,ee in m:
            if col in sidx: e[sidx[col]]+=ee
            else: t=t*pow(x0[col],ee,P)%P
        if t:
            key=tuple(e); out[key]=(out.get(key,0)+t)%P
    return {m:c for m,c in out.items() if c}
polys=[(i,sub(f)) for i,f in lin]+[(i,sub(f if not isinstance(f,int) else {():f})) for i,f in non]
varying=[(i,g) for i,g in polys if any(any(e) for e in g)]
aff=[(i,g) for i,g in varying if max(sum(m) for m in g)<=1]
hi=[(i,g) for i,g in varying if max(sum(m) for m in g)>1]
rowsA=[]
for i,g in aff:
    r=[0]*k
    for m,c in g.items():
        if sum(m)==1: r[[j for j,e in enumerate(m) if e][0]]=c%P
    rowsA.append((i,tuple(r)))
grp=collections.defaultdict(list)
for i,r in rowsA:
    j0=min(j for j in range(k) if r[j]%P); iv=pow(r[j0],-1,P)
    grp[tuple(x*iv%P for x in r)].append(i)
dirs=list(grp.items())
print('affine rows %d -> distinct directions %d ; multiplicities %s'
      %(len(rowsA),len(dirs),sorted(len(v) for _,v in dirs)))
def rank(rows):
    M=[list(r) for r in rows]; rr=0
    for c in range(k):
        pr=None
        for t in range(rr,len(M)):
            if M[t][c]%P: pr=t;break
        if pr is None: continue
        M[rr],M[pr]=M[pr],M[rr]
        iv=pow(M[rr][c],-1,P); M[rr]=[x*iv%P for x in M[rr]]
        for t in range(len(M)):
            if t!=rr and M[t][c]%P:
                f=M[t][c]; M[t]=[(x-f*y)%P for x,y in zip(M[t],M[rr])]
        rr+=1
    return rr
R_all=rank([d for d,_ in dirs])
print('rank of all directions: %d (unknowns %d, cost-free kernel %d)'%(R_all,k,k-R_all))
def kernel(rows):
    M=[list(r) for r in rows]; piv=[]; rr=0
    for c in range(k):
        pr=None
        for t in range(rr,len(M)):
            if M[t][c]%P: pr=t;break
        if pr is None: continue
        M[rr],M[pr]=M[pr],M[rr]
        iv=pow(M[rr][c],-1,P); M[rr]=[x*iv%P for x in M[rr]]
        for t in range(len(M)):
            if t!=rr and M[t][c]%P:
                f=M[t][c]; M[t]=[(x-f*y)%P for x,y in zip(M[t],M[rr])]
        piv.append(c); rr+=1
    fr=[c for c in range(k) if c not in piv]
    K=[]
    for c in fr:
        vec=[0]*k; vec[c]=1
        for t,cp in enumerate(piv): vec[cp]=(-M[t][c])%P
        K.append(vec)
    return K
base=[x0[c] for c in SUP]
ctx=flint.fmpz_mod_poly_ctx(P)
def cubic_ok(K):
    """can all higher-degree equations be zeroed on base + span(K)?  K includes the
    cost-free kernel; the cubics depend only on the extra directions."""
    m=len(K)
    def expand(g):
        out=collections.defaultdict(int)
        for mono,c in g.items():
            cur={(0,)*m:c%P}
            for j,e in enumerate(mono):
                for _ in range(e):
                    linj={(0,)*m:base[j]}
                    for t in range(m):
                        ee=[0]*m; ee[t]=1
                        if K[t][j]%P: linj[tuple(ee)]=K[t][j]%P
                    nxt=collections.defaultdict(int)
                    for m1,c1 in cur.items():
                        for m2,c2 in linj.items():
                            nxt[tuple(a+b for a,b in zip(m1,m2))]=(nxt[tuple(a+b for a,b in zip(m1,m2))]+c1*c2)%P
                    cur={a:b for a,b in nxt.items() if b}
            for a,b in cur.items(): out[a]=(out[a]+b)%P
        return {a:b for a,b in out.items() if b}
    E=[expand(g) for _,g in hi]
    nz=[e for e in E if e]
    if not nz: return True,0
    # count how many are nonzero constants (cannot be fixed)
    constbad=[e for e in nz if not any(any(a) for a in e)]
    if constbad: return False,len(constbad)
    # try: solve by univariate gcd if only one effective parameter
    eff=[t for t in range(m) if any(a[t] for e in nz for a in e)]
    if len(eff)==1:
        t=eff[0]; g=None
        for e in nz:
            deg=max(a[t] for a in e)
            f=ctx([sum(c for a,c in e.items() if a[t]==d)%P for d in range(deg+1)])
            g=f if g is None else g.gcd(f)
        return (g.degree()>0 and bool(g.roots())), 0
    if len(eff)==2:
        t1,t2=eff
        # bivariate: eliminate t2 by Sylvester resultant of the two lowest-degree polys,
        # then test every root of the resultant against all equations
        def as_t2(e):
            d2=max(a[t2] for a in e)
            return [ {a[t1]:c for a,c in e.items() if a[t2]==j} for j in range(d2+1) ]
        polys2=sorted(nz,key=lambda e:max(a[t2] for a in e))
        import itertools as _it
        cands=set()
        for e in polys2:
            d2=max(a[t2] for a in e)
            if d2==0:
                dd=max(a[t1] for a in e)
                f=ctx([sum(c for a,c in e.items() if a[t1]==j)%P for j in range(dd+1)])
                for r,_m in (f.roots() or []): cands.add(int(r))
        if not cands:
            # resultant of the first two in t2
            import sympy
            x,y=sympy.symbols('x y')
            def sp(e):
                return sum(int(c)*x**a[t1]*y**a[t2] for a,c in e.items())
            try:
                R=sympy.Poly(sympy.resultant(sp(polys2[0]),sp(polys2[1]),y),x)
                co=[int(c)%P for c in R.all_coeffs()][::-1]
                f=ctx(co)
                for r,_m in (f.roots() or []): cands.add(int(r))
            except Exception:
                return None,len(eff)
        for r1 in cands:
            g=None
            for e in nz:
                d2=max(a[t2] for a in e)
                co=[sum(c*pow(r1,a[t1],P) for a,c in e.items() if a[t2]==j)%P for j in range(d2+1)]
                while co and co[-1]==0: co.pop()
                if not co: continue
                if len(co)==1: g=ctx([1]); break
                f=ctx(co); g=f if g is None else g.gcd(f)
            if g is not None and g.degree()>0 and g.roots(): return True,0
        return False,len(eff)
    return None,len(eff)
t0=time.time(); found=[]
for w in range(1,BUDGET+1):
    for combo in itertools.combinations(range(len(dirs)),1 if w==1 else None) if False else []:
        pass
    # enumerate subsets of directions with total multiplicity exactly <= w
    idxs=list(range(len(dirs)))
    def rec(start,chosen,tot):
        if tot>w: return
        if chosen:
            rest=[dirs[j][0] for j in idxs if j not in chosen]
            if rank(rest)<R_all:
                K=kernel(rest)
                viol=sorted(i for j in chosen for i in dirs[j][1])
                ok,info=cubic_ok(K)
                print('   violated %d eqs %s -> cubics zeroable: %s (%s)'%(len(viol),viol,ok,info),flush=True)
                if ok: found.append((len(viol),viol))
        for j in range(start,len(idxs)):
            m=len(dirs[j][1])
            if tot+m<=w: rec(j+1,chosen+[j],tot+m)
    rec(0,[],0)
    if found:
        print('\n*** minimum total cost with all cubics zeroed: %d ***'%min(f[0] for f in found)); break
    print('  budget %d: no admissible violated-set (%.0fs)'%(w,time.time()-t0),flush=True)
