import heal_harness as H
from collections import defaultdict
p=H.p
d=H.loadd('best_agentA_39022.json')
base={v:d.get(v,0) for v in H.freeinp}
base[4287]=1; base[2081]=1
# keep x_9118, x_8731 at 0 baseline to measure affine dependence cleanly
base[9118]=0; base[8731]=0
for v in H.freeinp: H.val[v]=base[v]
H.forward(); V=H.val
desc_of=defaultdict(list)
for k,t in enumerate(H.order):
    for w in H.anc[t]: desc_of[w].append(k)
ns={'v':V,'__builtins__':{}}
touch9118=sorted(set(desc_of[9118])); touch8731=sorted(set(desc_of[8731]))
GQ=[2239,31731,9106,35619,23754,9629,18253]
def setf(v,x,tt):
    V[v]=x
    for k in tt: V[H.order[k]]=eval(H.gcode[k],ns)
# baseline values (x_8731=x_9118=0)
b={q:V[q] for q in GQ}
print("baseline (x_8731=x_9118=0):",{('x_%d'%q):str(b[q])[:30] for q in GQ})
# slope wrt x_9118
setf(9118,1,touch9118); s9118={q:V[q]-b[q] for q in GQ}; setf(9118,0,touch9118)
# slope wrt x_8731
setf(8731,1,touch8731); s8731={q:V[q]-b[q] for q in GQ}; setf(8731,0,touch8731)
# verify affine: check at (x_9118=3,x_8731=5)
setf(9118,3,touch9118); setf(8731,5,touch8731)
ok=all(V[q]==b[q]+3*s9118[q]+5*s8731[q] for q in GQ)
setf(9118,0,touch9118); setf(8731,0,touch8731)
print("affine verified?",ok)
for q in GQ:
    print(f"x_{q} = {b[q]} + ({s9118[q]})*x_9118 + ({s8731[q]})*x_8731")
# check which are wires
print("\nx_30095=%d (=p? %s)"%(V[30095],V[30095]==p))
print("x_24490=%d (=p? %s)"%(V[24490],V[24490]==p))
print("x_26874=%d (=p? %s)"%(V[26874],V[26874]==p))
import json
json.dump({'b':{q:str(b[q]) for q in GQ},'s9118':{q:str(s9118[q]) for q in GQ},'s8731':{q:str(s8731[q]) for q in GQ},
           'x7068':str(d.get(7068,0)),'x4432':str(d.get(4432,0))},
          open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/gadget_affine.json','w'))
