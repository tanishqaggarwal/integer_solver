#!/usr/bin/env python3
"""Support cone for ALL sensitive atoms (deg-1 + deg-2 + deg-4). Classify boolean constraints."""
import sg2_tl as T
import pickle
from collections import deque
p = T.p
sens = pickle.load(open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/sens.pkl','rb'))['sens']
# constraints = all sensitive atoms; targets G1,G2 (20862,20864); 42669 slaved -> treat as constraint too? it's currently nonzero.
constraints = list(sens)
for g in (20862,20864):
    if g not in constraints: constraints.append(g)
print(f"constraints (all sensitive): {len(constraints)}")

gate_inputs = {}
for t in T.order:
    _, rhs, vids = T.gates[T.definer[t]]
    gate_inputs[t] = tuple(vids)

checkvars = set()
for ai in constraints: checkvars |= T.L.atom_vars(T.A[ai]['poly'])
relevant=set(); support_free=set()
q=deque(checkvars); seen=set(checkvars)
while q:
    v=q.popleft()
    if v in T.freeinp: support_free.add(v); continue
    relevant.add(v)
    for u in gate_inputs.get(v,()):
        if u not in seen: seen.add(u); q.append(u)
rel_order=[t for t in T.order if t in relevant]
print(f"checkvars {len(checkvars)}, relevant gates {len(relevant)}, support free {len(support_free)}")
pickle.dump({'constraints':constraints,'relevant':relevant,'support_free':support_free,
             'rel_order':rel_order,'gate_inputs':gate_inputs},
            open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/support2.pkl','wb'))
