import os,sys,json,time,random
os.chdir('/home/user/integer_solver/solve_lab'); sys.path.insert(0,'.')
import heal_harness as H
from collections import defaultdict
p=H.p
pins=json.load(open('pinrec.json'))
allbits=set(r[1] for r in pins)
vA=H.loadd('best_agentA_39022.json')
base={v:vA.get(v,0) for v in H.freeinp}
for v in H.freeinp: H.val[v]=base[v]
H.val[13195]=1; H.forward()
desc_of=defaultdict(list)
for k,t in enumerate(H.order):
    for w in H.anc[t]: desc_of[w].append(k)
ns={'v':H.val,'__builtins__':{}}
def setf(v,x):
    H.val[v]=x
    for k in desc_of[v]: H.val[H.order[k]]=eval(H.gcode[k],ns)
def mfails(): return [i for i,c in enumerate(H.eqcode) if eval(c,ns)%p!=0]
# preset active-bit pin targets (free) to CONST
for atomidx,sel,target,const,coef,handle in pins:
    if H.val[sel]%p==1 and target in H.freeinp:
        setf(target, const%p)
print("after pin preset:",len(mfails()),"mod-p fails")
def inv(x): return pow(x%p,p-2,p)
def newton_step(Fp):
    K=set()
    for i in Fp:
        for v in H.eqvars[i]: K|=(H.anc.get(v,{v})&H.freeinp)
    K-=allbits; K=sorted(K); nK=len(K)
    b=[eval(H.eqcode[i],ns)%p for i in Fp]
    bc={v:H.val[v] for v in K}
    m=len(Fp)
    Aug=[[0]*(nK+1) for _ in range(m)]
    for j,v in enumerate(K):
        o=bc[v]; setf(v,o+1)
        for r in range(m): Aug[r][j]=(eval(H.eqcode[Fp[r]],ns)-b[r])%p
        setf(v,o)
    for r in range(m): Aug[r][nK]=(-b[r])%p
    # gaussian to get particular solution (free=0), track consistency
    pr=0; pivcol=[-1]*m
    for c in range(nK):
        sel=None
        for r in range(pr,m):
            if Aug[r][c]%p: sel=r;break
        if sel is None: continue
        Aug[pr],Aug[sel]=Aug[sel],Aug[pr]
        iv=inv(Aug[pr][c]); Aug[pr]=[(x*iv)%p for x in Aug[pr]]
        for r in range(m):
            if r!=pr and Aug[r][c]%p:
                f=Aug[r][c]; Aug[r]=[(Aug[r][t]-f*Aug[pr][t])%p for t in range(nK+1)]
        pivcol[pr]=c; pr+=1
        if pr==m: break
    sol=[0]*nK
    for r in range(pr):
        sol[pivcol[r]]=Aug[r][nK]%p
    # apply
    for j,v in enumerate(K):
        if sol[j]: setf(v,(bc[v]+sol[j])%p)
    # residual reduction count
    return
t0=time.time(); hist=[]
for it in range(40):
    Fp=mfails(); hist.append(len(Fp))
    if not Fp: break
    newton_step(Fp)
    if it<8 or it%5==0: print(f"iter {it}: fails={len(Fp)} t={time.time()-t0:.0f}s",flush=True)
    if time.time()-t0>1500: print("timebound"); break
print("history:",hist)
if hist and hist[-1]==0:
    print("*** 0 mod-p fails! saving ***")
    json.dump({('x_%d'%v):str(H.val[v]) for v in H.freeinp}, open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/modp13195b.json','w'))
