#!/usr/bin/env python3
"""Parse product monomials (x_i*x_j) in the constraint equations. Identify which of the moved
handles form 'conflict pairs' (both factors of a product both move) causing 2nd-order breakage.
Goal: find product-active free inputs so we can pick CLEAN compensators."""
import json, ast, re, time
from agentC_common import (p, gates, order, definer, forward, val, freeinp, ns, lines, eqcode,
                           eqvars, load_best, CORE, downstream_ks, partial_forward, rootcode_of, inv)
from collections import defaultdict

best=load_best(); forward()
V=json.load(open('agentC_Vdata.json'))
H=V['H']; Hset=set(H)

# --- parse product monomials in equations: collect pairs of variable-ids multiplied together ---
# We consider the gate DEFINITIONS too, since equations reference gates that expand to products.
# For breakage, what matters is degree-2 monomials over FREE INPUTS after full expansion. That's
# expensive; instead detect product atoms syntactically in eq text and gate rhs, then map factors
# to their free-input cones (a product breaks 2nd-order iff both factors' cones contain moved frees).
gate_defs={t:(rhs,vids) for t,rhs,vids in gates}
def freecone(root, memo={}):
    if root in memo: return memo[root]
    seen=set(); leaves=set(); st=[root]
    while st:
        x=st.pop()
        if x in seen: continue
        seen.add(x)
        if x in gate_defs:
            for u in gate_defs[x][1]: st.append(u)
        elif x in freeinp: leaves.add(x)
    memo[root]=leaves; return leaves

VAR=re.compile(r'x_(\d+)')
def product_pairs(expr):
    """Return list of (setA,setB) variable-id sets that are multiplied (each side's vars)."""
    node=ast.parse(expr,mode='eval').body
    pairs=[]
    def vids(n): return set(int(m) for m in VAR.findall(ast.unparse(n)))
    def walk(n):
        if isinstance(n,ast.BinOp):
            if isinstance(n.op,ast.Mult):
                la=vids(n.left); lb=vids(n.right)
                # only care if both sides contain variables (not constant*var)
                if la and lb:
                    pairs.append((la,lb))
            walk(n.left); walk(n.right)
        elif isinstance(n,ast.UnaryOp):
            walk(n.operand)
    walk(node)
    return pairs

# the 14 broken eqs from root0
broken=[5814, 8388, 12678, 13589, 13813, 14682, 16372, 17904, 22492, 23342, 28207, 33744, 34857, 37277]
# moved handles for root0: all H with nonzero null combo. Compute delta.
Hidx={h:i for i,h in enumerate(H)}; nullbasis=V['nullbasis']
x29322=val[29322]%p; x1326=val[1326]%p; x27713=val[27713]%p; x33469=val[33469]%p; x3558=val[3558]%p
T0=val[6671]%p; den=(x29322+x1326)%p
beta=((x27713+x3558)*inv(den))%p; alpha=((-T0)*inv(den))%p
import agentC_poly as Ply
f1=Ply.pmul([x33469,1],Ply.pmul([x29322,(-1)%p],[x29322,(-1)%p]))
lin=[(x3558-alpha)%p,(-beta)%p]; Spoly=Ply.psub(f1,Ply.pmul(lin,lin))
roots=Ply.roots_mod_p(Spoly)
da=roots[0]; dc=(alpha+beta*da)%p
delta={h:(da*nullbasis[0][Hidx[h]]+dc*nullbasis[1][Hidx[h]])%p for h in H}
moved=set(h for h in H if delta[h]!=0)
print(f"moved handles (root0): {len(moved)} of {len(H)}")

# For each broken eq, find product atoms where BOTH factor-cones contain moved frees
conflict_frees=set()
for i in broken:
    expr=lines[i].rsplit('=',1)[0]
    pairs=product_pairs(expr)
    for (la,lb) in pairs:
        ca=set(); cb=set()
        for v in la: ca|=freecone(v)
        for v in lb: cb|=freecone(v)
        ma=ca&moved; mb=cb&moved
        if ma and mb:
            conflict_frees|=ma|mb
print(f"conflict free inputs (appear as moved factors in broken products): {len(conflict_frees)}")
print(sorted(conflict_frees)[:40])
# are controls among conflicts?
print("controls in conflict:", [c for c in (12186,14853,16742) if c in conflict_frees])

# Which moved handles are 'product-active' anywhere in cons (appear in a product with a var whose
# cone has other moved frees)? Classify all H by product activity to pick clean compensators.
# Build global product-partner info over ALL equations is expensive; restrict to cons eqs.
eqbyvar=defaultdict(set)
for i in range(len(lines)):
    for v in eqvars[i]: eqbyvar[v].add(i)
F0=set(i for i in range(len(lines)) if eval(eqcode[i],ns)!=0)
cons=set()
for h in H:
    aff=set(eqbyvar.get(h,()))
    for k in downstream_ks(h): aff|=eqbyvar.get(order[k],set())
    cons|=(aff-F0)
print(f"|cons|={len(cons)}; scanning products in cons...", flush=True)
# free input -> set of product-partner free inputs (via cones) within cons
partner=defaultdict(set)
t0=time.time()
for ci,i in enumerate(cons):
    expr=lines[i].rsplit('=',1)[0]
    for (la,lb) in product_pairs(expr):
        ca=set(); cb=set()
        for v in la: ca|=freecone(v)
        for v in lb: cb|=freecone(v)
        # restrict to handles
        cah=ca&Hset; cbh=cb&Hset
        for x in cah:
            partner[x]|=cbh-{x}
        for x in cbh:
            partner[x]|=cah-{x}
print(f"scanned in {time.time()-t0:.0f}s")
prod_active=set(h for h in H if partner[h])
print(f"product-active handles (have a handle product-partner in cons): {len(prod_active)}")
print(f"product-inactive (clean) handles: {len(H)-len(prod_active)}")
print("controls product-active?", {c:(c in prod_active) for c in (12186,14853,16742)})
# partners of the controls
for c in (12186,16742):
    print(f"  x_{c} product-partners (handles): {sorted(partner[c])[:20]}")
json.dump({'prod_active':sorted(prod_active),
           'clean':sorted(set(H)-prod_active),
           'partner':{str(k):sorted(v) for k,v in partner.items() if v}},
          open('agentC_products.json','w'))
print("saved agentC_products.json")
