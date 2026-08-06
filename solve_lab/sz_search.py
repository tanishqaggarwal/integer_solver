import heal_harness as H, sz_engine as E, sz_inner as SI
import re,time,random,json
from collections import defaultdict
p=H.p; RIP=list(E.RIP); RIPS=set(E.RIP)
lines=[L for L in open('../EQUATIONS.txt').read().split('\n') if L.strip()]
PRODFAC={29638:1,35935:1}
E.classify(); ns={'v':H.val,'__builtins__':{}}
def inner_code(i):
    lhs=lines[i].rsplit('=',1)[0]; facs=SI.toplevel_factors(lhs); vf=[f for f in facs if '_' in f]
    base=vf[PRODFAC[i]] if i in PRODFAC else vf[0]
    return compile(re.sub(r'x_(\d+)',r'v[\1]',base),'<i>','eval')
IC={i:inner_code(i) for i in RIP}

# cone knobs
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

t=time.time(); 
d=H.loadd('best_agentA_39022.json')
for v in H.freeinp: H.val[v]=d.get(v,0)
H.forward(); f=time.time()-t
t=time.time(); nf=len(H.fails()); tf=time.time()-t
print(f"39022 baseline fails={nf}; forward+fails timing: fails()={tf*1000:.0f}ms")

def total_fails(): return len(H.fails())
def setup_hybrid():
    for v in H.freeinp: H.val[v]=d.get(v,0)
    H.forward(); r7=H.val[2099]; r4=H.val[19964]
    H.val[7068]=r7; H.val[4432]=r4; H.forward()

BEST=[nf, None]   # (fails, snapshot of changed frees)
def consider(tag):
    F=total_fails()
    if F<BEST[0]:
        BEST[0]=F; BEST[1]={v:H.val[v] for v in H.freeinp if H.val[v]!=d.get(v,0)}
        print(f"  *** NEW BEST {F} ({tag})")
    return F

# ---- Greedy exact single-eq-zero descent from hybrid ----
def greedy(maxpass=6):
    improved=True; pas=0
    while improved and pas<maxpass:
        improved=False; pas+=1
        F0=total_fails()
        curfail=[i for i in RIP if eval(IC[i],ns)!=0]
        for i in curfail:
            Ei=eval(IC[i],ns)
            if Ei==0: continue
            # find a knob with coef dividing Ei; try zeroing
            random.shuffle(knobs)
            for w in knobs[:40]:
                H.val[w]+=1
                for _ in (0,): H.forward()
                Ei1=eval(IC[i],ns); a=Ei1-Ei
                H.val[w]-=1; H.forward()
                if a!=0 and Ei%a==0:
                    step=-Ei//a
                    old=H.val[w]; H.val[w]=old+step; H.forward()
                    if eval(IC[i],ns)==0 and total_fails()<=F0:
                        nf=consider(f"greedy eq{i} knob{w}")
                        if nf<F0: improved=True; F0=nf; break
                        # keep if not worse
                    else:
                        H.val[w]=old; H.forward()
        # end pass
    return total_fails()

random.seed(1)
print("\n=== greedy from hybrid ===")
setup_hybrid(); print("hybrid fails",total_fails())
g=greedy(); print("after greedy:",g)

print("\n=== multi-start: random large perturb of a few cone knobs, then greedy ===")
results=[]
for s in range(12):
    setup_hybrid()
    for w in random.sample(knobs, min(4,len(knobs))):
        H.val[w]+= random.randrange(-p,p)
    H.forward()
    r=greedy(maxpass=3)
    results.append(r)
    print(f" start {s}: reached {r} fails (best so far {BEST[0]})")

# ---- also try from 39022 side: activate x_17325,x_9413 handles to shift G1/G2, measure ----
print("\n=== 39022 side: can coset moves reduce the 11? ===")
for v in H.freeinp: H.val[v]=d.get(v,0)
H.forward(); base_fail=total_fails()
# x_17325 shifts x_642 -> G1 in steps of 7376877*p ; sweep small k
best_local=base_fail
for k in range(-3,4):
    H.val[17325]=k; H.forward(); F=total_fails()
    if F<best_local: best_local=F
H.val[17325]=0
for k in range(-3,4):
    H.val[9413]=k; H.forward(); F=total_fails()
    if F<best_local: best_local=F
H.val[9413]=0; H.forward()
print(f"39022 with x_17325/x_9413 sweeps: best={best_local} (baseline {base_fail})")

print(f"\n==== GLOBAL BEST fails = {BEST[0]} (current known best 11) ====")
if BEST[0]<11 and BEST[1] is not None:
    out=dict(d); 
    for k,v in BEST[1].items(): out['x_%d'%k]=v
    # write full state
    full={('x_%d'%i):H.val[i] for i in range(H.NVARS)}
    json.dump({('x_%d'%k):int(v) for k,v in {**{int(kk[2:]):vv for kk,vv in d.items()}, **BEST[1]}.items()}, open('sz_best.json','w'))
    print("checkpointed sz_best.json")
else:
    print("no improvement over 11; nothing checkpointed")
