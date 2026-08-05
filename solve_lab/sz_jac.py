import heal_harness as H, sz_engine as E, sz_inner as SI
import re,time,sys
from collections import defaultdict
p=H.p; RIP=E.RIP; CORE=set(E.CORE); G1G2=set(E.G1G2)
lines=[L for L in open('../EQUATIONS.txt').read().split('\n') if L.strip()]
PRODFAC={29638:1,35935:1}   # zero-factor index for product ripple eqs
def inner_code(i):
    lhs=lines[i].rsplit('=',1)[0]
    facs=SI.toplevel_factors(lhs); vf=[f for f in facs if '_' in f]
    if i in PRODFAC: base=vf[PRODFAC[i]]
    elif len(vf)==2 and vf[0]==vf[1]: base=vf[0]
    elif len(vf)==1: base=vf[0]
    else:
        # generic prod: pick factor that is 0 at compensator (reference)
        base=vf[-1]
    return compile(re.sub(r'x_(\d+)',r'v[\1]',base),'<i>','eval')

E.classify(); r7,r4=E.setup()
ns={'v':H.val,'__builtins__':{}}
assert set(H.fails())==set(RIP)

# cone + knobs
cone=set(); stack=[]
for e in RIP:
    for w in H.eqvars[e]:
        if w not in cone: cone.add(w); stack.append(w)
while stack:
    w=stack.pop(); gi=H.definer.get(w)
    if gi is None: continue
    for u in H.gates[gi][2]:
        if u not in cone: cone.add(u); stack.append(u)
knobs=sorted(w for w in cone if w in H.freeinp)
print(f"knobs (cone frees): {len(knobs)}")

# blast radius T: eqs touched by knobs
setk=set(knobs)
desc=defaultdict(set)
for f in knobs: desc[f].add(f)
for t in H.order:
    for f in (H.anc[t]&setk): desc[f].add(t)
var_eqs=defaultdict(list)
for i,vs in enumerate(H.eqvars):
    for v in vs: var_eqs[v].append(i)
T=set()
for f in knobs:
    for t in desc[f]: T.update(var_eqs.get(t,()))
T=sorted(T)
print(f"blast radius |T|={len(T)}  (incl 16 ripple). core in T? {bool(set(T)&CORE)} g1g2 in T? {bool(set(T)&G1G2)}")

# inner codes + baseline residuals
ic={i:inner_code(i) for i in T}
def inner_vec():
    return {i:eval(ic[i],ns) for i in T}
r0=inner_vec()
nfail0=sum(1 for i in T if r0[i]!=0)
print(f"currently nonzero inner in T: {nfail0} (should be 16)")

# Jacobian mod p
t0=time.time()
col={w:j for j,w in enumerate(knobs)}; nc=len(knobs)
rows_map=defaultdict(dict)  # eq -> {col: coef}
base=H.val[:]
for w in knobs:
    H.val[w]+=1; H.forward()
    rv={i:eval(ic[i],ns) for i in T}
    for i in T:
        dcoef=(rv[i]-r0[i])%p
        if dcoef: rows_map[i][col[w]]=dcoef
    H.val[:]=base[:]
H.forward()
print(f"Jacobian built {time.time()-t0:.1f}s; rows with content: {sum(1 for i in T if rows_map[i])}")

# Build augmented system: for each eq in T, sum coef*delta = -r0[i] mod p
A=[]; rhs=[]; ids=[]
for i in T:
    row=[0]*nc
    for c,v in rows_map[i].items(): row[c]=v
    A.append(row); rhs.append((-r0[i])%p); ids.append(i)
# Gaussian elim mod p
def inv(a): return pow(a%p,p-2,p)
pr=0; where=[-1]*nc; m=len(A)
for c in range(nc):
    sel=-1
    for r in range(pr,m):
        if A[r][c]%p: sel=r; break
    if sel<0: continue
    A[pr],A[sel]=A[sel],A[pr]; rhs[pr],rhs[sel]=rhs[sel],rhs[pr]; ids[pr],ids[sel]=ids[sel],ids[pr]
    iv=inv(A[pr][c]); A[pr]=[(x*iv)%p for x in A[pr]]; rhs[pr]=(rhs[pr]*iv)%p
    for r in range(m):
        if r!=pr and A[r][c]%p:
            f=A[r][c]
            A[r]=[(A[r][k]-f*A[pr][k])%p for k in range(nc)]
            rhs[r]=(rhs[r]-f*rhs[pr])%p
    where[c]=pr; pr+=1
incons=[ids[r] for r in range(m) if all(A[r][k]%p==0 for k in range(nc)) and rhs[r]%p!=0]
print(f"rank={pr}  #inconsistent rows={len(incons)}")
if incons:
    inR=[i for i in incons if i in set(RIP)]; inO=[i for i in incons if i not in set(RIP)]
    print(f"  inconsistent among 16-ripple: {len(inR)} {inR}")
    print(f"  inconsistent among OTHER(satisfied): {len(inO)} sample {inO[:12]}")
    print("VERDICT: first-order heal over cone knobs is mod-p INCONSISTENT")
else:
    print("VERDICT: first-order CONSISTENT mod p over cone knobs!")
