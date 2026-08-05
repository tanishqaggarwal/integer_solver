import os,sys,json,time
os.chdir('/home/user/integer_solver/solve_lab'); sys.path.insert(0,'.')
import heal_harness as H
import flint
from collections import defaultdict
p=H.p; ctx=flint.fmpz_mod_ctx(p)
vA=H.loadd('best_agentA_39022.json')
base={v:vA.get(v,0) for v in H.freeinp}
for v in H.freeinp: H.val[v]=base[v]
H.val[13195]=1  # the forced bit
H.forward()
desc_of=defaultdict(list)
for k,t in enumerate(H.order):
    for w in H.anc[t]: desc_of[w].append(k)
var_eqs=defaultdict(list)
for i,vs in enumerate(H.eqvars):
    for v in vs: var_eqs[v].append(i)
ns={'v':H.val,'__builtins__':{}}
def setf(v,x):
    H.val[v]=x
    for k in desc_of[v]: H.val[H.order[k]]=eval(H.gcode[k],ns)
def mfails():
    return [i for i,c in enumerate(H.eqcode) if eval(c,ns)%p!=0]
def rref_solve(rows_entries, nfp, K):
    # rows_entries: list over fp of dict{col:val}, plus rhs in col nK
    nK=len(K)
    M=flint.fmpz_mod_mat(nfp, nK+1, ctx)
    for r,(ent,rhs) in enumerate(rows_entries):
        for c,val in ent.items(): M[r,c]=val
        M[r,nK]=rhs%p
    Mr=M.rref()[0] if isinstance(M.rref(),tuple) else M.rref()
    # find pivot per row, build particular solution (free=0)
    sol=[0]*nK
    incons=False
    for r in range(nfp):
        piv=None
        for c in range(nK):
            if int(Mr[r,c])%p!=0: piv=c;break
        if piv is None:
            if int(Mr[r,nK])%p!=0: incons=True
        else:
            sol[piv]=int(Mr[r,nK])  # since rref, pivot=1, rhs is the value with free=0
    return sol,incons
t0=time.time()
hist=[]
for it in range(25):
    Fp=mfails()
    hist.append(len(Fp))
    if not Fp: break
    # free ancestors of Fp
    K=set()
    for i in Fp:
        for v in H.eqvars[i]: K|= (H.anc.get(v,{v}) & H.freeinp)
    K.discard(13195)
    K=sorted(K)
    kidx={v:j for j,v in enumerate(K)}
    # residual b
    b={i:eval(H.eqcode[i],ns)%p for i in Fp}
    # jacobian rows
    rows=[]
    bc={v:H.val[v] for v in K}
    # build per-col deltas
    coldelta=defaultdict(dict)
    for v in K:
        o=bc[v]; setf(v,o+1)
        for i in Fp:
            d=(eval(H.eqcode[i],ns)-b[i])%p
            if d: coldelta[i][kidx[v]]=d
        setf(v,o)
    rows_entries=[(coldelta[i], (-b[i])%p) for i in Fp]
    sol,incons=rref_solve(rows_entries,len(Fp),K)
    if incons:
        print(f"iter {it}: {len(Fp)} fails, linear system INCONSISTENT mod p")
        break
    for j,v in enumerate(K):
        if sol[j]: setf(v,(bc[v]+sol[j])%p)
    print(f"iter {it}: fails={len(Fp)}, |K|={len(K)}, t={time.time()-t0:.0f}s")
print("history:",hist)
if hist and hist[-1]==0:
    print("CONVERGED to 0 mod-p fails!")
    json.dump({('x_%d'%v):str(H.val[v]%p) for v in H.freeinp}, open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/modp_sol.json','w'))
