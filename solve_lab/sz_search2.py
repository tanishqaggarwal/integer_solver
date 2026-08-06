import heal_harness as H, sz_engine as E, sz_inner as SI
import re,time,random,json
from math import gcd
p=H.p; RIP=list(E.RIP)
lines=[L for L in open('../EQUATIONS.txt').read().split('\n') if L.strip()]
PRODFAC={29638:1,35935:1}
E.classify(); ns={'v':H.val,'__builtins__':{}}
def inner_code(i):
    lhs=lines[i].rsplit('=',1)[0]; facs=SI.toplevel_factors(lhs); vf=[f for f in facs if '_' in f]
    base=vf[PRODFAC[i]] if i in PRODFAC else vf[0]
    return compile(re.sub(r'x_(\d+)',r'v[\1]',base),'<i>','eval')
IC={i:inner_code(i) for i in RIP}
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
d=H.loadd('best_agentA_39022.json')
def hybrid():
    for v in H.freeinp: H.val[v]=d.get(v,0)
    H.forward(); H.val[7068]=H.val[2099]; H.val[4432]=H.val[19964]; H.forward()
def fails(): return len(H.fails())
def egcdx(a,b):
    old_r,r=a,b; old_s,s=1,0; old_t,t=0,1
    while r!=0:
        q=old_r//r
        old_r,r=r,old_r-q*r
        old_s,s=s,old_s-q*s
        old_t,t=t,old_t-q*t
    return (old_r,old_s,old_t)
def jac_row(i,ks):
    Ei=eval(IC[i],ns); out={}
    for w in ks:
        H.val[w]+=1; H.forward(); out[w]=eval(IC[i],ns)-Ei; H.val[w]-=1
    H.forward(); return Ei,out
BEST=[11,None]
def consider(tag):
    F=fails()
    if F<BEST[0]:
        BEST[0]=F; BEST[1]={v:H.val[v] for v in H.freeinp if H.val[v]!=d.get(v,0)}
        print(f"   *** NEW BEST {F} ({tag})")
    return F
def gs_round(ks, tries=16):
    curfail=[i for i in RIP if eval(IC[i],ns)!=0]
    random.shuffle(curfail)
    for i in curfail[:tries]:
        Ei,row=jac_row(i,ks)
        if Ei==0: continue
        cand=[(w,c) for w,c in row.items() if c!=0]
        if len(cand)<2: continue
        (w1,a1),(w2,a2)=cand[0],cand[1]
        g,x0,y0=egcdx(abs(a1),abs(a2))
        if Ei% g!=0: 
            # try other pairs for gcd|Ei
            done=False
            for j in range(len(cand)):
                for k in range(j+1,min(j+4,len(cand))):
                    w1,a1=cand[j]; w2,a2=cand[k]; g,x0,y0=egcdx(abs(a1),abs(a2))
                    if Ei%g==0: done=True;break
                if done:break
            if not done: continue
        s1=1 if a1>0 else -1; s2=1 if a2>0 else -1
        m=-Ei//g
        d1=s1*x0*m; d2=s2*y0*m
        o1,o2=H.val[w1],H.val[w2]; F0=fails()
        H.val[w1]=o1+d1; H.val[w2]=o2+d2; H.forward()
        if fails()<=F0 and eval(IC[i],ns)==0:
            consider(f"gs eq{i}")
        else:
            H.val[w1]=o1; H.val[w2]=o2; H.forward()
random.seed(7)
T0=time.time()
print("knobs",len(knobs))
hybrid(); print("hybrid fails",fails())
for r in range(9):
    if time.time()-T0>210: print("time cap"); break
    hybrid()
    if r>0:
        for w in random.sample(knobs,min(3,len(knobs))): H.val[w]+=random.randrange(-p,p)
        H.forward()
    for _ in range(3): gs_round(random.sample(knobs,min(24,len(knobs))))
    print(f" restart {r}: fails now {fails()}  BEST {BEST[0]}  ({time.time()-T0:.0f}s)")
print(f"\n==== SEARCH BEST = {BEST[0]} (known best 11) ====")
if BEST[0]<11 and BEST[1]:
    full={int(k[2:]) if str(k).startswith('x_') else int(k):v for k,v in d.items()}
    full.update(BEST[1])
    json.dump({('x_%d'%k):int(v) for k,v in full.items()}, open('sz_best.json','w'))
    print("checkpointed sz_best.json with",BEST[0],"fails")
else:
    print("no improvement; the wall holds at 11")
