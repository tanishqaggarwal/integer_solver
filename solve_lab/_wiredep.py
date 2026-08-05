import heal_harness as H
p=H.p
# which free vars does each wire depend on?
wires=[17499,28599,26874,22665,28961,7075,13859,15616]
d=H.loadd('best_agentA_39022.json')
for v in H.freeinp: H.val[v]=d.get(v,0)
H.forward(); V=H.val
knobs=[7068,4432,17325,9413,2099,19964]
for w in wires:
    if w not in H.anc:
        print(f"x_{w}: (free input)"); continue
    fa=H.anc[w]; nf=len(fa)
    onknobs = sorted(fa & set(knobs))
    print(f"x_{w}=... (=p? {V[w]==p}): #free_anc={nf}, depends on knobs={onknobs}")
# Why is x_9413 nonlinear in residuals? Check: does x_9413 feed a wire that multiplies a free var?
# Check ancestors overlap: is x_9413 an ancestor of x_17499?
print("\nx_9413 in anc(x_17499)?", 9413 in H.anc.get(17499,set()))
print("x_17325 in anc(x_28599)?", 17325 in H.anc.get(28599,set()))
# The nonlinearity source: perturb x_9413, does x_17499 change?
from collections import defaultdict
desc_of=defaultdict(list)
for k,t in enumerate(H.order):
    for wv in H.anc[t]: desc_of[wv].append(k)
ns={'v':V,'__builtins__':{}}
def setfree(vv,x):
    V[vv]=x
    for k in desc_of[vv]: V[H.order[k]]=eval(H.gcode[k],ns)
b17499=V[17499]; b28599=V[28599]; b28730=V[28730]; b642=V[642]
setfree(9413, d.get(9413,0)+1)
print(f"\nafter x_9413+=1: x_17499 changed? {V[17499]!=b17499} (delta bits {abs(V[17499]-b17499).bit_length()}); x_28730 delta = {V[28730]-b28730} (=p? {V[28730]-b28730==p})")
setfree(9413, d.get(9413,0))
setfree(17325, d.get(17325,0)+1)
print(f"after x_17325+=1: x_28599 changed? {V[28599]!=b28599}; x_642 delta = {V[642]-b642} (=p? {V[642]-b642==p})")
setfree(17325, d.get(17325,0))
