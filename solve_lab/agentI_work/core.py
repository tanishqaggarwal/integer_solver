"""Neighborhood of the 7 residual atoms: all atoms touching their variables."""
import sys, collections, os
from model import Model, load_assign
M = Model()
v = load_assign(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                             '..', 'best', 'new_instance_partial_39026.json'))
RES = [23432, 23433, 36225, 36226, 36227, 36228, 36229]
var2atoms = collections.defaultdict(list)
for i, vs in enumerate(M.avars):
    for x in vs:
        var2atoms[x].append(i)

seed = set()
for a in RES:
    seed |= M.avars[a]
print("seed vars:", sorted(seed))
print()
for x in sorted(seed):
    print(f"--- X{x} = {v[x]!r}"[:200])
    for a in var2atoms[x]:
        val = M.atom_val(a, v)
        s = f"    a{a}: {M.src[a]}"
        if len(s) > 160: s = s[:160] + "..."
        print(s, " val=", ("0" if val == 0 else "NONZERO"), f" eqs={len(M.atom_eqs[a])}")
