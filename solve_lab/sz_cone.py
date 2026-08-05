import heal_harness as H, sz_engine as E
import json,time
from collections import defaultdict
p=H.p; RIP=set(E.RIP); CORE=set(E.CORE); G1G2=set(E.G1G2)
E.classify(); r7,r4=E.setup()
assert set(H.fails())==RIP

# ripple cone
cone=set(); stack=[]
for e in E.RIP:
    for w in H.eqvars[e]:
        if w not in cone: cone.add(w); stack.append(w)
while stack:
    w=stack.pop(); gi=H.definer.get(w)
    if gi is None: continue
    for u in H.gates[gi][2]:
        if u not in cone: cone.add(u); stack.append(u)
cone_free=sorted(w for w in cone if w in H.freeinp)
print(f"cone {len(cone)}  cone_free {len(cone_free)}")

# global fanout: eq -> free deps (only need reverse for cone_free). Build var->eqs, then free->eqs via descendants.
t0=time.time()
# descendants of each cone_free among gate vars: t s.t. cone_free-member in anc[t]
setcf=set(cone_free)
desc=defaultdict(set)   # free -> set of gate/leaf vars it feeds (incl itself)
for f in cone_free: desc[f].add(f)
for t in H.order:
    inter=H.anc[t]&setcf
    for f in inter: desc[f].add(t)
# var -> eqs
var_eqs=defaultdict(list)
for i,vs in enumerate(H.eqvars):
    for v in vs: var_eqs[v].append(i)
def fanout_eqs(f):
    s=set()
    for t in desc[f]: s.update(var_eqs.get(t,()))
    return s
print(f"built desc in {time.time()-t0:.1f}s")

# classify each cone free by its fanout
private=[]      # touches only RIP
ripplus=[]      # touches RIP + others (all others currently satisfied? we know only RIP fail)
for f in cone_free:
    fo=fanout_eqs(f)
    touches_core = bool(fo&CORE)
    touches_g = bool(fo&G1G2)
    others = fo - RIP
    if not others:
        private.append(f)
    else:
        ripplus.append((f,len(fo),len(others),touches_core,touches_g))
print(f"\nPRIVATE cone frees (touch ONLY the 16 ripple): {len(private)}")
print(private)
print(f"\nnon-private cone frees: {len(ripplus)}  (f, nfanout, n_other_eqs, hitsCORE, hitsG1G2)")
for row in sorted(ripplus,key=lambda x:x[2])[:25]: print(row)
# how many non-private touch NEITHER core nor g1g2 (only ripple + 'inert' satisfied eqs)
safe_extra=[r for r in ripplus if not r[3] and not r[4]]
print(f"\nnon-private but avoid core&G1G2: {len(safe_extra)}")
json.dump({'private':private,'cone_free':cone_free},open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/sz_cone.json','w'))
