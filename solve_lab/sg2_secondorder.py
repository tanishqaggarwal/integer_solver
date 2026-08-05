#!/usr/bin/env python3
"""Search for second-order DOF: product monomials (a,b) in sensitive atoms where BOTH factors
are currently 0 (invisible to first-order Jacobian). Also products wire*free where free!=0.
These are candidate 'missing DOF' that could shift the invariant K."""
import sg2_lib as L, pickle
import heal_harness as H
p=H.p
atoms=L.load_atoms_full(); A={a['idx']:a for a in atoms}
sens=pickle.load(open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/sens.pkl','rb'))['sens']
gaps=pickle.load(open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/gaps.pkl','rb'))
wire=gaps['wire']
vA=H.loadd('best_agentA_39022.json')
for v in H.freeinp: H.val[v]=vA.get(v,0)
H.forward()
val=H.val

# gate map: which product gates exist
prod_gate={}
for t in H.order:
    _,rhs,vids=H.gates[H.definer[t]]
    if rhs.count('*')==1 and '+' not in rhs and '-' not in rhs:
        parts=[int(x[2:]) for x in rhs.replace(' ','').split('*') if x.startswith('x_')]
        if len(parts)==2: prod_gate[t]=tuple(parts)

# Among sensitive atoms, find deg-2 monomials (a,b) with val[a]==val[b]==0
sensset=set(sens)
double_zero=[]  # (atom, a, b)
for ai in sens:
    for m,c in A[ai]['poly'].items():
        if len(m)==2:
            a,b=m
            if val[a]==0 and val[b]==0:
                double_zero.append((ai,a,b,c))
print(f"double-zero product monomials in sensitive atoms: {len(double_zero)}")
# classify: how many involve a free input among the pair
freeinp=H.freeinp
dz_free=[(ai,a,b) for ai,a,b,c in double_zero if a in freeinp or b in freeinp]
print(f"  ... with at least one FREE factor: {len(dz_free)}")
# unique factor vars that are free and zero
zfree=set()
for ai,a,b,c in double_zero:
    if a in freeinp and val[a]==0: zfree.add(a)
    if b in freeinp and val[b]==0: zfree.add(b)
print(f"  distinct free zero factors: {len(zfree)}")

# Now look specifically at the THREE gap atoms' product structure and their upstream relay
# The gaps: 20862(G1),20864(G2),41390(G3). Their p-slacks: x_642,x_28730,x_24410 (wire*free).
# Check: do any UPSTREAM relay gaps (the 71 in the conservation law) have NON-wire product slacks?
# Load the 71-combo from cert? recompute quickly: just scan all 604 gap atoms for non-wire product terms.
gapatoms=[g[0] for g in gaps['gaps']]
nonwire_prod_gaps=[]
for aidx in gapatoms:
    for m,c in A[aidx]['poly'].items():
        if len(m)==1:
            t=m[0]
            if t in prod_gate:
                a,b=prod_gate[t]
                if a not in wire and b not in wire:
                    nonwire_prod_gaps.append((aidx,t,a,b))
print(f"gap atoms with a NON-wire product slack term: {len(nonwire_prod_gaps)}")
for x in nonwire_prod_gaps[:20]:
    aidx,t,a,b=x
    print(f"  atom {aidx}: pslack x_{t}=x_{a}(v={val[a]%p!=0})*x_{b}(v={val[b]%p!=0})  {A[aidx]['repr'][:60]}")
