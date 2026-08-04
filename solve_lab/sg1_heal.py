"""Set residue-fixed + handles state, then damped coordinate descent using ALL frees as
compensators (Gauss-Seidel), quadratic-root candidates, monotone fail-count acceptance."""
import sys, json, time, math, random
sys.path.insert(0,'/home/user/integer_solver/solve_lab')
import heal_harness as H
p=H.p; val=H.val; ns={'v':val,'__builtins__':{}}
eqcode=H.eqcode; eqvars=H.eqvars
from collections import defaultdict
efs=[]
for i in range(len(eqvars)):
    s=set()
    for var in eqvars[i]:
        if var in H.freeinp: s.add(var)
        else: s|=H.anc.get(var,set())
    efs.append(s)
free_to_eqs=defaultdict(set)
for i,s in enumerate(efs):
    for f in s: free_to_eqs[f].add(i)

def setup():
    d=H.loadd('best/new_instance_partial_39013.json')
    for v in H.freeinp: val[v]=d.get(v,0)
    H.forward()
    r29=val[29322]%p; r35=val[3558]%p
    val[14853]-=r29; val[16742]+=r35
    H.forward()
    L1=val[11150]; L2=val[25739]; L3=val[37758]
    val[30317]=(-L1)//p; val[2936]=537773*L3//p; val[5146]=L2//(6672769*p)
    H.forward()
setup()
allfails=lambda: set(i for i,c in enumerate(eqcode) if eval(c,ns)!=0)
F=allfails()
print(f"residue-fixed+handles: {len(F)} fails")

def candidates_for(f,relF):
    cands=set(); x0=val[f]; samp={}
    for dx in (0,1,2):
        val[f]=x0+dx; H.forward()
        samp[dx]=[eval(eqcode[i],ns) for i in relF]
    val[f]=x0; H.forward()
    for k,i in enumerate(relF):
        y0,y1,y2=samp[0][k],samp[1][k],samp[2][k]
        d2=y2-2*y1+y0
        if d2%2!=0:
            if y1!=y0 and y0%(y1-y0)==0: cands.add(x0-y0//(y1-y0))
            continue
        A=d2//2; B=y1-y0-A; C=y0
        if A==0:
            if B!=0 and (-C)%B==0: cands.add(x0+(-C)//B)
        else:
            disc=B*B-4*A*C
            if disc>=0:
                r=math.isqrt(disc)
                if r*r==disc:
                    for s in (r,-r):
                        num=-B+s
                        if num%(2*A)==0: cands.add(x0+num//(2*A))
    return cands

def total_res(Fs): return sum(abs(eval(eqcode[i],ns)) for i in Fs)

# knobs: all frees in fail support (checked+unchecked)
def greedy(F,budget=180):
    t0=time.time(); best=len(F)
    while time.time()-t0<budget:
        knobs=set()
        for i in F: knobs|=efs[i]
        knobs=list(knobs); random.shuffle(knobs); improved=False
        for f in knobs:
            relF=list(F&free_to_eqs[f])
            if not relF: continue
            cands=candidates_for(f,relF)
            x0=val[f]; bx=x0; bn=len(F); br=total_res(F); bF=F
            for x in cands:
                if x==x0: continue
                val[f]=x; H.forward()
                nF=set(F)
                for i in free_to_eqs[f]:
                    if eval(eqcode[i],ns)!=0: nF.add(i)
                    else: nF.discard(i)
                n=len(nF)
                if n<bn or (n==bn and total_res(nF)<br):
                    bn=n; bx=x; br=total_res(nF); bF=nF
            val[f]=bx; H.forward()
            if bx!=x0 and bn<len(F):
                F=bF; improved=True
                if len(F)<best:
                    best=len(F); print(f"  new best {best} (f={f}) t={time.time()-t0:.0f}s",flush=True)
        if not improved:
            print(f"  plateau {len(F)}",flush=True); break
    return F

F=greedy(F)
print(f"FINAL: {len(F)} fails: {sorted(F)}")
if len(F)<=22:
    json.dump({f"x_{i}":val[i] for i in range(H.NVARS)}, open('sg1_heal_out.json','w'))
    print("saved sg1_heal_out.json")
