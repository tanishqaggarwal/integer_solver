import os,sys,json,time
os.chdir('/home/user/integer_solver/solve_lab'); sys.path.insert(0,'.')
import heal_harness as H
from collections import defaultdict
p=H.p
vA=H.loadd('best_agentA_39022.json')
base={v:vA.get(v,0) for v in H.freeinp}
for v in H.freeinp: H.val[v]=base[v]
H.forward()
FAILS11=[2554,6816,8124,8680,9421,12231,12270,12350,14584,22044,29125]
RIPPLE16=[697,1985,5225,10815,16048,17784,17801,22402,23667,24721,27124,28737,29638,29959,35935,37431]
LOCAL=set(FAILS11+RIPPLE16)
bits10=[2081,4287,5910,11368,13195,17406,18022,22562,23751,28005]
Eanc=set()
for i in LOCAL:
    for v in H.eqvars[i]: Eanc|=H.anc.get(v,{v})
cont=sorted(v for v in Eanc if v in H.freeinp and v not in bits10)
desc_of=defaultdict(list)
for k,t in enumerate(H.order):
    for w in H.anc[t]: desc_of[w].append(k)
var_eqs=defaultdict(list)
for i,vs in enumerate(H.eqvars):
    for v in vs: var_eqs[v].append(i)
ns={'v':H.val,'__builtins__':{}}
def setf_full(v,x):
    H.val[v]=x
    for k in desc_of[v]: H.val[H.order[k]]=eval(H.gcode[k],ns)
# for each cont input, which NON-local eqs break when perturbed?
F0=set(H.fails())
safe=[]; unsafe=[]
for v in cont:
    o=base[v]
    te=set(var_eqs[v])
    for k in desc_of[v]: te.update(var_eqs[H.order[k]])
    b0={i:eval(H.eqcode[i],ns) for i in te}
    setf_full(v,o+1)
    broke=set(i for i in te if (eval(H.eqcode[i],ns)==0)!=(b0[i]==0))
    setf_full(v,o)
    nonlocal_break = broke - LOCAL
    if nonlocal_break: unsafe.append((v,len(nonlocal_break)))
    else: safe.append(v)
print(f"cont free inputs: {len(cont)}")
print(f"SAFE (only affect local 27 eqs): {len(safe)}")
print(f"  {safe}")
print(f"UNSAFE (also break non-local eqs): {len(unsafe)}")
for v,n in sorted(unsafe,key=lambda x:x[1])[:20]:
    print(f"  x_{v}: breaks {n} non-local eqs")
json.dump({'safe':safe,'cont':cont,'unsafe':[u[0] for u in unsafe]}, open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/safe.json','w'))
