import os,sys,json,time
os.chdir('/home/user/integer_solver/solve_lab'); sys.path.insert(0,'.')
import heal_harness as H
import flint
from collections import defaultdict
p=H.p; ctx=flint.fmpz_mod_ctx(p)
pins=json.load(open('pinrec.json'))
allbits=set(r[1] for r in pins)
vA=H.loadd('best_agentA_39022.json')
base={v:vA.get(v,0) for v in H.freeinp}
for v in H.freeinp: H.val[v]=base[v]
H.val[13195]=1
H.forward()
desc_of=defaultdict(list)
for k,t in enumerate(H.order):
    for w in H.anc[t]: desc_of[w].append(k)
ns={'v':H.val,'__builtins__':{}}
def setf(v,x):
    H.val[v]=x
    for k in desc_of[v]: H.val[H.order[k]]=eval(H.gcode[k],ns)
def mfails():
    return [i for i,c in enumerate(H.eqcode) if eval(c,ns)%p!=0]
def rref_particular(Fp,K):
    kidx={v:j for j,v in enumerate(K)}; nK=len(K)
    b={i:eval(H.eqcode[i],ns)%p for i in Fp}
    coldelta=defaultdict(dict); bc={v:H.val[v] for v in K}
    for v in K:
        o=bc[v]; setf(v,o+1)
        for i in Fp:
            d=(eval(H.eqcode[i],ns)-b[i])%p
            if d: coldelta[i][kidx[v]]=d
        setf(v,o)
    M=flint.fmpz_mod_mat(len(Fp),nK+1,ctx)
    for r,i in enumerate(Fp):
        for c,val in coldelta[i].items(): M[r,c]=val
        M[r,nK]=(-b[i])%p
    R=M.rref()
    R=R[0] if isinstance(R,tuple) else R
    sol=[0]*nK; incons=False
    for r in range(len(Fp)):
        piv=None
        for c in range(nK):
            if int(R[r,c])%p!=0: piv=c;break
        if piv is None:
            if int(R[r,nK])%p!=0: incons=True
        else: sol[piv]=int(R[r,nK])%p
    return sol,incons,bc
t0=time.time(); hist=[]
for it in range(30):
    Fp=mfails(); hist.append(len(Fp))
    if not Fp: break
    K=set()
    for i in Fp:
        for v in H.eqvars[i]: K|=(H.anc.get(v,{v})&H.freeinp)
    K-=allbits   # EXCLUDE all boolean bits
    K=sorted(K)
    sol,incons,bc=rref_particular(Fp,K)
    if incons:
        print(f"iter {it}: {len(Fp)} fails, |K|={len(K)}, INCONSISTENT")
        break
    for j,v in enumerate(K):
        if sol[j]: setf(v,(bc[v]+sol[j])%p)
    if it<6 or it%5==0: print(f"iter {it}: fails={len(Fp)} |K|={len(K)} t={time.time()-t0:.0f}s")
print("history:",hist)
if hist and hist[-1]==0:
    print("*** CONVERGED to 0 mod-p fails! ***")
    json.dump({('x_%d'%v):str(H.val[v]) for v in H.freeinp}, open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/modp13195.json','w'))
