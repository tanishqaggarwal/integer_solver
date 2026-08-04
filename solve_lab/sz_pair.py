import heal_harness as H, sz_engine as E, sz_inner as SI
import re,time
from collections import defaultdict
p=H.p; RIP=set(E.RIP); CORE=set(E.CORE); G1G2=set(E.G1G2)
lines=[L for L in open('../EQUATIONS.txt').read().split('\n') if L.strip()]
PRODFAC={29638:1,35935:1}
def inner_code(i):
    lhs=lines[i].rsplit('=',1)[0]; facs=SI.toplevel_factors(lhs); vf=[f for f in facs if '_' in f]
    base = vf[PRODFAC[i]] if i in PRODFAC else (vf[0] if len(vf)>=1 else lhs)
    return compile(re.sub(r'x_(\d+)',r'v[\1]',base),'<i>','eval')
E.classify(); E.setup(); ns={'v':H.val,'__builtins__':{}}
assert set(H.fails())==RIP

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
# topo index + children map (parents = definer vids)
topo={t:k for k,t in enumerate(H.order)}
children=defaultdict(list)
for t in H.order:
    for v in H.gates[H.definer[t]][2]:
        children[v].append(t)
# strict descendants (defined gate vars) of a cone var, topo-sorted
def desc_topo(c):
    seen=set(); st=[c]; out=[]
    while st:
        x=st.pop()
        for t in children.get(x,()):
            if t not in seen: seen.add(t); st.append(t)
    return sorted(seen,key=lambda t:topo[t])
# blast radius T (eqs touched by cone) via var_eqs on cone descendants
var_eqs=defaultdict(list)
for i,vs in enumerate(H.eqvars):
    for v in vs: var_eqs[v].append(i)
Tset=set()
for c in cone: Tset.update(var_eqs.get(c,()))
T=sorted(Tset)
print(f"cone {len(cone)}  |T| {len(T)}")
ic={i:inner_code(i) for i in T}
r0={i:eval(ic[i],ns) for i in T}
print("nonzero inner in T:",sum(1 for i in T if r0[i]!=0))
# INJECT-any-cone-var Jacobian (surgical). knobs = all cone vars (gates+frees).
knobs=[c for c in cone if H.definer.get(c) is not None or c in H.freeinp]  # all cone vars
dcache={c:desc_topo(c) for c in knobs}
base=H.val[:]
col={c:j for j,c in enumerate(knobs)}; nc=len(knobs)
rows=defaultdict(dict)
t0=time.time()
for c in knobs:
    H.val[c]+=1
    for t in dcache[c]: H.val[t]=eval(H.gcode[topo[t]],ns)
    for i in T:
        dv=(eval(ic[i],ns)-r0[i])%p
        if dv: rows[i][col[c]]=dv
    H.val[:]=base[:]
print(f"inject-Jacobian built {time.time()-t0:.1f}s over {nc} cone vars")
# consistency: zero the 16, keep others
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
inR=[i for i in incons if i in RIP]; inO=[i for i in incons if i not in RIP]
print(f"MAX-INJECT rank={pr}  inconsistent={len(incons)}  ripple-blocked={len(inR)} {sorted(inR)}  other-blocked={len(inO)}")
if not inR:
    print(">>> The 16 ripple ARE reachable via cone-internal injection (realizability then matters)")
else:
    print(">>> Even arbitrary cone-internal injection CANNOT zero these ripple: rigid obstruction")
