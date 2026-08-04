import heal_harness as H
import json
from collections import defaultdict
p=H.p
vA=H.loadd('best_agentA_39022.json')
for v in H.freeinp: H.val[v]=vA.get(v,0)
H.forward()
V=H.val[:]
baseF=set(H.fails())
print("baseline fails:", sorted(baseF))
ns={'v':H.val,'__builtins__':{}}
def eq_resid(i): return eval(H.eqcode[i],ns)
# For a set of free-input perturbations, forward and report which SATISFIED eqs break
def probe(perturb, label):
    for v in range(len(H.val)): H.val[v]=V[v]
    for v,dv in perturb.items(): H.val[v]+=dv
    H.forward()
    newF=set(H.fails())
    broke=sorted(newF-baseF)
    fixed=sorted(baseF-newF)
    print(f"{label}: total={len(newF)} broke={len(broke)} fixed={len(fixed)}")
    if broke: print("   broke:",broke[:30])
    # restore
    for v in range(len(H.val)): H.val[v]=V[v]
    H.forward()
    return newF

# Which eqs does each free knob appear in, and are they satisfied?
knobs=[7068,4432,17325,9413,6947,26874,23754,14853,12186,16742]
for k in knobs:
    if k not in H.freeinp:
        print(f"x_{k}: NOT free"); continue
    eqs=[i for i,vs in enumerate(H.eqvars) if k in vs]
    sat=[i for i in eqs if i not in baseF]
    print(f"x_{k}: in {len(eqs)} eqs, {len(sat)} satisfied (would break): {sat[:15]}")
