import heal_harness as H
import ast
from collections import defaultdict
p=H.p
d=H.loadd('best/new_instance_partial_39013.json')
for v in H.freeinp: H.val[v]=d.get(v,0)
H.forward()

# Parse each gate rhs (string) into a poly: dict monomial(tuple sorted) -> coef
# gates in H order: gcode compiled; but we need the source rhs string. Reload gates.jsonl.
import json
gates={}
with open('atoms/gates.jsonl') as f:
    for line in f:
        dd=json.loads(line)
        gates.setdefault(dd['t'],[]).append((dd['rhs'],tuple(dd['vids'])))
# H.order gives the chosen target order; the definer rhs is what forward used.
# Reconstruct definer rhs from heal_harness: it recompiled gcode from gates[definer[t]].
# Simpler: parse the rhs string used. heal_harness stored gcode but not source. Re-derive:
# heal_harness picks first ready gate per target via its topo. We'll parse ALL candidate rhs and
# trust that forward() matches; verify by re-eval.
import re
VAR=re.compile(r'x_(\d+)')
def parse_poly(expr):
    node=ast.parse(expr,mode='eval').body
    def go(n):
        if isinstance(n,ast.Constant): return {():n.value}
        if isinstance(n,ast.Name): return {(int(n.id[2:]),):1}
        if isinstance(n,ast.UnaryOp): 
            return {m:-c for m,c in go(n.operand).items()}
        a=go(n.left); b=go(n.right)
        if isinstance(n.op,ast.Add):
            o=defaultdict(int)
            for m,c in a.items(): o[m]+=c
            for m,c in b.items(): o[m]+=c
            return {m:c for m,c in o.items() if c}
        if isinstance(n.op,ast.Sub):
            o=defaultdict(int)
            for m,c in a.items(): o[m]+=c
            for m,c in b.items(): o[m]-=c
            return {m:c for m,c in o.items() if c}
        if isinstance(n.op,ast.Mult):
            o=defaultdict(int)
            for m1,c1 in a.items():
                for m2,c2 in b.items():
                    o[tuple(sorted(m1+m2))]+=c1*c2
            return {m:c for m,c in o.items() if c}
        raise ValueError
    return go(node)

# Build gate poly per target in H.order using the source from heal_harness's definer.
# heal_harness: definer[t]=gi index into its own 'gates' list (same file). Reconstruct identically.
gsrc=[]
with open('atoms/gates.jsonl') as f:
    for line in f:
        dd=json.loads(line); gsrc.append((dd['t'],dd['rhs'],tuple(dd['vids'])))
# replicate topo from heal_harness to get definer rhs per order target
# Actually H.gcode[k] is compiled from gates[definer[order[k]]][1]; we can get the source by
# matching: but easier—just eval poly at base and compare to H.val to validate parse.
# We'll parse using H's order + a re-run of the same selection. Instead: for each t in H.order,
# find a candidate rhs whose poly evaluates (at current val) to H.val[t]. Use that.
def polyeval(poly):
    s=0
    for m,c in poly.items():
        tt=c%p
        for v in m: tt=(tt*H.val[v])%p
        s=(s+tt)%p
    return s
cand=defaultdict(list)
for (t,rhs,vids) in gsrc: cand[t].append((rhs,vids))
gatepoly={}
mismatch=0
for t in H.order:
    chosen=None
    for rhs,vids in cand[t]:
        pl=parse_poly(rhs)
        if polyeval(pl)==H.val[t]%p:
            chosen=pl; break
    if chosen is None:
        chosen=parse_poly(cand[t][0][0]); mismatch+=1
    gatepoly[t]=chosen
print("gate polys built. mismatches:", mismatch, "of", len(H.order))
# degrees
from collections import Counter
degc=Counter()
for t,pl in gatepoly.items():
    degc[max((len(m) for m in pl),default=0)]+=1
print("gate max-degree histogram:", dict(degc))
# G = sum over gates of |free ancestors|
G=sum(len(H.anc[t]) for t in H.order)
print("G = sum_gates |free anc| =", G)
# sum over gates of |inputs|
GI=sum(len(H.anc[t]) for t in H.order)  # placeholder
import pickle
pickle.dump(gatepoly, open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/gatepoly.pkl','wb'))
