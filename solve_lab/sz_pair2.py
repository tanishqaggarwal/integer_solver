import heal_harness as H, sz_engine as E, sz_inner as SI
import re,time
from collections import defaultdict
p=H.p; RIP=set(E.RIP)
lines=[L for L in open('../EQUATIONS.txt').read().split('\n') if L.strip()]
PRODFAC={29638:1,35935:1}
E.classify(); E.setup(); ns={'v':H.val,'__builtins__':{}}
assert set(H.fails())==RIP

def factor_codes(i):
    lhs=lines[i].rsplit('=',1)[0]; facs=SI.toplevel_factors(lhs); vf=[f for f in facs if '_' in f]
    return [compile(re.sub(r'x_(\d+)',r'v[\1]',f),'<i>','eval') for f in vf]
def choose_inner(i):
    fc=factor_codes(i)
    if i in RIP:
        if i in PRODFAC: return fc[PRODFAC[i]]
        return fc[0]                       # square base / lin form
    # satisfied eq: pick the factor currently == 0 (min |val|)
    vals=[abs(eval(c,ns)) for c in fc]
    return fc[vals.index(min(vals))]

# cone
cone=set(); stack=[]
for e in RIP:
    for w in H.eqvars[e]:
        if w not in cone: cone.add(w); stack.append(w)
while stack:
    w=stack.pop(); gi=H.definer.get(w)
    if gi is None: continue
    for u in H.gates[gi][2]:
        if u not in cone: cone.add(u); stack.append(u)
topo={t:k for k,t in enumerate(H.order)}
children=defaultdict(list)
for t in H.order:
    for v in H.gates[H.definer[t]][2]: children[v].append(t)
def desc_topo(c):
    seen=set(); st=[c]
    while st:
        x=st.pop()
        for t in children.get(x,()):
            if t not in seen: seen.add(t); st.append(t)
    return sorted(seen,key=lambda t:topo[t])
var_eqs=defaultdict(list)
for i,vs in enumerate(H.eqvars):
    for v in vs: var_eqs[v].append(i)
Tset=set()
for c in cone: Tset.update(var_eqs.get(c,()))
T=sorted(Tset)
ic={i:choose_inner(i) for i in T}
r0={i:eval(ic[i],ns) for i in T}
nz=[i for i in T if r0[i]!=0]
print(f"cone {len(cone)} |T| {len(T)}  nonzero-inner {len(nz)}  (==16 ripple? {set(nz)==RIP})")
if set(nz)!=RIP: print("  extra nonzero:",[i for i in nz if i not in RIP][:10])

knobs=list(cone)
dcache={c:desc_topo(c) for c in knobs}
base=H.val[:]; col={c:j for j,c in enumerate(knobs)}; nc=len(knobs)
rows=defaultdict(dict)
t0=time.time()
for c in knobs:
    H.val[c]+=1
    for t in dcache[c]: H.val[t]=eval(H.gcode[topo[t]],ns)
    for i in T:
        dv=(eval(ic[i],ns)-r0[i])%p
        if dv: rows[i][col[c]]=dv
    H.val[:]=base[:]
print(f"inject-Jac {time.time()-t0:.1f}s; constrained eqs (rows w content or r0!=0): building")
A=[]; rhs=[]; ids=[]
for i in T:
    if not rows[i] and r0[i]%p==0: continue
    row=[0]*nc
    for cc,v in rows[i].items(): row[cc]=v
    A.append(row); rhs.append((-r0[i])%p); ids.append(i)
def inv(a): return pow(a%p,p-2,p)
pr=0; m=len(A)
for c in range(nc):
    sel=-1
    for r in range(pr,m):
        if A[r][c]%p: sel=r;break
    if sel<0: continue
    A[pr],A[sel]=A[sel],A[pr]; rhs[pr],rhs[sel]=rhs[sel],rhs[pr]; ids[pr],ids[sel]=ids[sel],ids[pr]
    iv=inv(A[pr][c]); A[pr]=[(x*iv)%p for x in A[pr]]; rhs[pr]=(rhs[pr]*iv)%p
    for r in range(m):
        if r!=pr and A[r][c]%p:
            f=A[r][c]; A[r]=[(A[r][k]-f*A[pr][k])%p for k in range(nc)]; rhs[r]=(rhs[r]-f*rhs[pr])%p
    pr+=1
incons=[ids[r] for r in range(m) if all(A[r][k]%p==0 for k in range(nc)) and rhs[r]%p!=0]
inR=sorted(i for i in incons if i in RIP); inO=[i for i in incons if i not in RIP]
print(f"MAX-INJECT (correct factors): rank={pr} m={m}  ripple-blocked={len(inR)} {inR}  other-satisfied-blocked={len(inO)} {inO[:12]}")
print(">>> RIGID" if inR else ">>> ripple reachable by cone injection; realizability is the remaining question")
