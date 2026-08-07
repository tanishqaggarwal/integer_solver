"""Search for a departure direction with affine cost <= 6 on which the 20 cubic
equations can be zeroed.  Randomised: a weight-w codeword is orthogonal to 65-w rows,
so a random rank-10 subset of the rows lies in its kernel with probability ~0.36."""
import os, sys, pickle, random, itertools, collections, time
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentG_work')
import gsym2 as G
from gsym2 import L, ad, P
import flint
D=pickle.load(open('/home/user/integer_solver/solve_lab/agentG_work/coset_model.pkl','rb'))
Lin=pickle.load(open('/home/user/integer_solver/solve_lab/agentG_work/coset_lin.pkl','rb'))
NB=D['NB']; lin=D['lin']; non=D['non']; x0=Lin['x0']; pt=D['pt']
ix={u:i for i,u in enumerate(NB)}
SUP=[ix[int(u)] for u in sys.argv[1].split(',')]
TRIALS=int(sys.argv[2]) if len(sys.argv)>2 else 3000
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
    rowsA.append((i,r))
E=len(rowsA); base=[x0[c] for c in SUP]
def kernel(rows):
    M=[r[:] for r in rows]; piv=[]; rr=0
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
    return K,rr
K0,r0=kernel([r for _,r in rowsA])
print('all rows: rank %d, cost-free kernel dim %d'%(r0,len(K0)))
def cost(v): return sum(1 for _,r in rowsA if sum(a*b for a,b in zip(r,v))%P)
# cubics on the line base + s*v (the cost-free kernel leaves them constant)
def line_polys(v):
    out=[]
    for i,g in hi:
        co=collections.defaultdict(int)
        for mono,c in g.items():
            # expand prod_j (base_j + s v_j)^{e_j}
            cur={0:c%P}
            for j,e in enumerate(mono):
                for _ in range(e):
                    nxt=collections.defaultdict(int)
                    for dg,cc in cur.items():
                        if base[j]%P: nxt[dg]=(nxt[dg]+cc*base[j])%P
                        if v[j]%P: nxt[dg+1]=(nxt[dg+1]+cc*v[j])%P
                    cur={a:b for a,b in nxt.items() if b}
            for a,b in cur.items(): co[a]=(co[a]+b)%P
        out.append((i,{a:b for a,b in co.items() if b}))
    return out
ctx=flint.fmpz_mod_poly_ctx(P)
def common_root(lp):
    g=None
    for i,co in lp:
        if not co: continue
        deg=max(co); f=ctx([co.get(t,0) for t in range(deg+1)])
        g=f if g is None else g.gcd(f)
        if g.degree()==0: return None
    if g is None: return 'ANY'
    rts=g.roots()
    return rts if rts else None
d_star=[(pt[c]-x0[c])%P for c in SUP]
print('deliverable direction: affine cost %d ; cubics common root %s'%(cost(d_star),str(common_root(line_polys(d_star)))[:60]))
random.seed(1)
allrows=list(range(E))
best=(10**9,None); seen={}
t0=time.time()
for trial in range(TRIALS):
    S=random.sample(allrows,10)
    Ks,rr=kernel([rowsA[t][1] for t in S])
    if rr<10: continue
    # extra directions beyond the cost-free kernel
    for v in Ks:
        c=cost(v)
        if c==0: continue
        key=tuple(sorted((j,val) for j,val in enumerate(v) if val))
        if c<=8 and key not in seen:
            seen[key]=c
            lp=line_polys(v); rt=common_root(lp)
            viol=[i for i,r in rowsA if sum(a*b for a,b in zip(r,v))%P]
            print('  cost %d dir found (violates %s) cubics-common-root=%s'%(c,viol,str(rt)[:50]),flush=True)
            if rt and rt!='ANY' and c<best[0]: best=(c,(v,rt,viol))
    if trial%500==0: print('  trial %d, %d distinct low-cost directions, %.0fs'%(trial,len(seen),time.time()-t0),flush=True)
print('\ndistinct directions with affine cost <= 8 found: %d'%len(seen))
print('cost histogram:',dict(sorted(collections.Counter(seen.values()).items())))
print('BEST with a cubic common root: cost %s'%(best[0] if best[1] else 'none'))
