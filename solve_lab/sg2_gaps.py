#!/usr/bin/env python3
"""Enumerate ALL gap/relay atoms: find wire members (=p), product gates wire*free_rare (p-slacks),
and the atoms using them. Characterize the reduced gap system."""
import sg2_lib as L, pickle
import heal_harness as H
p=H.p
atoms=L.load_atoms_full(); A={a['idx']:a for a in atoms}
idx=pickle.load(open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/idx.pkl','rb'))
var_atoms=idx['var_atoms']
vA=H.loadd('best_agentA_39022.json')
for v in H.freeinp: H.val[v]=vA.get(v,0)
H.forward()

# wire members: gate vars with value == p (mod: exactly p)
wire=set(v for v in range(L.NVARS) if H.val[v]==p)
print(f"vars valued exactly p (wire members + others): {len(wire)}")

# product gates x_t = x_a * x_b
prod_gates={}
for t in H.order:
    _,rhs,vids=H.gates[H.definer[t]]
    if '*' in rhs and rhs.count('*')==1 and '+' not in rhs and '-' not in rhs:
        # simple product a*b
        parts=[int(x[2:]) for x in rhs.replace(' ','').split('*') if x.startswith('x_')]
        if len(parts)==2:
            prod_gates[t]=tuple(parts)

# p-slacks: product gate where one factor is wire (=p) and the other is a FREE input
pslacks={}  # t -> (partner_free, wire_factor)
for t,(a,b) in prod_gates.items():
    for f,w in ((a,b),(b,a)):
        if w in wire and f in H.freeinp:
            pslacks[t]=(f,w)
print(f"p-slack product gates (wire*free): {len(pslacks)}")
# how rare are the partners
partner_atoms={f:len(var_atoms[f]) for t,(f,w) in pslacks.items()}
import collections
print("partner #atoms distribution:", dict(collections.Counter(partner_atoms.values())))
# gaps: atoms (deg-1) containing a p-slack term with a FREE checked var
gaps=[]
for t,(f,w) in pslacks.items():
    for aidx in var_atoms[t]:
        a=A[aidx]
        if aidx in (H.definer.get(t) and []): pass
        deg=max((len(m) for m in a['poly']),default=0)
        # deg-1 atom that is a residue check (contains t linearly and a free var)
        if deg==1 and t in L.atom_vars(a['poly']):
            frees=[v for v in L.atom_vars(a['poly']) if v in H.freeinp and v!=f]
            if frees:
                gaps.append((aidx,t,f,w,frees))
# dedup by atom
seen=set(); ug=[]
for g in gaps:
    if g[0] in seen: continue
    seen.add(g[0]); ug.append(g)
print(f"gap atoms (deg-1 residue checks with p-slack + free checked var): {len(ug)}")
for aidx,t,f,w,frees in ug[:40]:
    val_res=H.val[frees[0]] if frees else 0
    print(f"  atom {aidx}: {A[aidx]['repr'][:75]}")
    print(f"      pslack x_{t}=x_{f}(free,{partner_atoms[t]}at)*x_{w}(wire); checked free {frees}")
pickle.dump({'wire':wire,'pslacks':pslacks,'gaps':ug}, open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/gaps.pkl','wb'))
