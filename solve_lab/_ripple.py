import heal_harness as H
from collections import defaultdict
p=H.p
d=H.loadd('best_agentA_39022.json')
for v in H.freeinp: H.val[v]=d.get(v,0)
H.forward()
F0=set(H.fails())
# equations each free var appears in (directly via eqvars)
eqof=defaultdict(list)
for i,vs in enumerate(H.eqvars):
    for v in vs: eqof[v].append(i)
targets=[6418,7068,4432,12553,17325,9413,2099,19964,642,28730]
for v in targets:
    eqs=eqof[v]
    infail=[i for i in eqs if i in F0]
    print(f"x_{v}: appears in {len(eqs)} equations directly; {len(infail)} are in the 11-fail set: {infail}")
# But eqvars only counts DIRECT appearance. The free vars also feed equations via gates.
# Find equations whose value depends on each free var (via descendants reaching eq vars).
desc=defaultdict(set)  # free var -> set of all variables it influences (itself + gate descendants)
for t in H.order:
    fa=H.anc[t]&H.freeinp
    for f in fa: desc[f].add(t)
for f in H.freeinp: desc[f].add(f)
# eq depends on free f if any eqvar is in desc[f]
eqdep=defaultdict(set)
var_to_eqs=defaultdict(list)
for i,vs in enumerate(H.eqvars):
    for v in vs: var_to_eqs[v].append(i)
print("\n--- full (through-gate) equation dependence ---")
for v in targets:
    influenced=desc[v]
    eqs=set()
    for w in influenced:
        eqs.update(var_to_eqs.get(w,[]))
    infail=sorted(eqs&F0); nonfail=sorted(eqs-F0)
    print(f"x_{v}: influences {len(eqs)} eqs total; {len(infail)} fail {infail}; {len(nonfail)} currently-SAT")
    if len(nonfail)<=40: print(f"    SAT eqs (ripple): {nonfail}")
