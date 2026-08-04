#!/usr/bin/env python3
"""Compute the backward cone (relevant gates, support free inputs) of the 606 constraint checks."""
import sg2_tl as T
import pickle
from collections import deque
p = T.p
sens = pickle.load(open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/sens.pkl','rb'))['sens']
deg1 = [ai for ai in sens if max((len(m) for m in T.A[ai]['poly']),default=0)==1]
constraints = deg1[:]  # 604 deg-1
# ensure G1(20862), G2(20864) included as targets
for g in (20862, 20864):
    if g not in constraints: constraints.append(g)
print(f"constraints: {len(constraints)} (incl G1,G2)")

# gate parents: for each gate target, its input vids
gate_inputs = {}
for t in T.order:
    _, rhs, vids = T.gates[T.definer[t]]
    gate_inputs[t] = tuple(vids)

# backward cone of all check vars
checkvars = set()
for ai in constraints:
    checkvars |= T.L.atom_vars(T.A[ai]['poly'])
relevant = set()  # gate vars in backward cone
support_free = set()
q = deque(checkvars)
seen = set(checkvars)
while q:
    v = q.popleft()
    if v in T.freeinp:
        support_free.add(v); continue
    relevant.add(v)
    for u in gate_inputs.get(v, ()):
        if u not in seen:
            seen.add(u); q.append(u)
print(f"checkvars: {len(checkvars)}")
print(f"relevant gate vars (backward cone): {len(relevant)}")
print(f"support free inputs: {len(support_free)}")
# topological order of relevant gates (subset of T.order)
rel_order = [t for t in T.order if t in relevant]
print(f"relevant gates in topo order: {len(rel_order)}")
pickle.dump({'constraints':constraints,'relevant':relevant,'support_free':support_free,
             'rel_order':rel_order,'gate_inputs':gate_inputs},
            open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/support.pkl','wb'))
