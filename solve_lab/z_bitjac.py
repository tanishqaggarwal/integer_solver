import os,sys,json,time
os.chdir('/home/user/integer_solver/solve_lab'); sys.path.insert(0,'.')
import heal_harness as H
import flint
from collections import defaultdict
p=H.p
ctx=flint.fmpz_mod_ctx(p)
pins=json.load(open('pinrec.json'))
selectors=sorted(set(r[1] for r in pins))
nb=len(selectors); bidx={b:j for j,b in enumerate(selectors)}
vA=H.loadd('best_agentA_39022.json')
base={v:vA.get(v,0) for v in H.freeinp}
for v in H.freeinp: H.val[v]=base[v]
H.forward()
desc_of=defaultdict(list)
for k,t in enumerate(H.order):
    for w in H.anc[t]: desc_of[w].append(k)
ns={'v':H.val,'__builtins__':{}}
var_eqs=defaultdict(list)
for i,vs in enumerate(H.eqvars):
    for v in vs: var_eqs[v].append(i)
def setfree(v,x):
    H.val[v]=x
    for k in desc_of[v]: H.val[H.order[k]]=eval(H.gcode[k],ns)
def touched_eqs(v):
    s=set(var_eqs[v])
    for k in desc_of[v]: s.update(var_eqs[H.order[k]])
    return s
# baseline residuals (only need touched eqs, but compute lazily)
t0=time.time()
# Build sparse Jacobian rows: for each bit, delta over its touched eqs
entries=defaultdict(dict)  # eq -> {bitcol: val}
allrows=set()
for b in selectors:
    old=base[b]; new=1-old if old in (0,1) else 0
    te=touched_eqs(b)
    b0={i:eval(H.eqcode[i],ns)%p for i in te}
    setfree(b,new)
    for i in te:
        d=(eval(H.eqcode[i],ns)-b0[i])%p
        if d!=0:
            entries[i][bidx[b]]=d; allrows.add(i)
    setfree(b,old)
print(f"jacobian built in {time.time()-t0:.0f}s; nonzero rows={len(allrows)}")
# build flint matrix rows x 256
rows=sorted(allrows)
M=flint.fmpz_mod_mat(len(rows), nb, ctx)
for ri,i in enumerate(rows):
    for c,val in entries[i].items():
        M[ri,c]=val
r=M.rank()
print(f"bit-Jacobian: {len(rows)} eqs x {nb} bits, rank={r}")
print(f"=> effective FREE bits (nullity) = {nb - r} = {nb} - {r}")
