"""S10 step 31: the canonical FORWARD-EVAL frame.

Take the 39,026 witness's FREE INPUT values, forward-evaluate every gate in
topological order (iterating for the ~1800 cyclic vars), and report which CHECK
atoms fail.  That is the honest global statement of the problem:
  choose 7,273 free inputs so that all 10,792 checks vanish.
"""
import os, sys, collections, json
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L
import tools as T

P = 2**256 - 2**32 - 977
base = L.load(os.path.join(LAB, 'best', 'new_instance_partial_39026.json'))

# --- what constrains the two free inputs that sit in the binding checks? ---
for u in (14853, 24548):
    print(f'\n=== x_{u} (FREE) occurs in {len(L.var_atoms[u])} atoms ===')
    for a in sorted(L.var_atoms[u]):
        out = L.atom_out.get(a)
        kind = f'GATE->x_{out[1]}' if out else 'CHECK'
        print(f'   a{a:<6} {kind:<14} neq={len(L.atom2eq.get(a,{})):<3} {L.atom_src[a][:110]}')

# --- forward evaluation --------------------------------------------------
print('\n=== forward evaluation from the witness free inputs ===')
topo = L.topo
definer, atom_out = L.definer, L.atom_out
v = list(base)
# reset every defined variable, then recompute
order = list(topo)
inorder = set(order)
rest = [x for x in definer if x not in inorder]
print(f'topo covers {len(order)}, cyclic remainder {len(rest)}')

def recompute(v, seq):
    for u in seq:
        a = definer[u]
        c, t = atom_out[a]
        nv = T.solve_lin(a, u, v)
        if nv is not None:
            v[u] = nv

for it in range(6):
    recompute(v, order)
    recompute(v, rest)
    av = L.all_atom_values(v)
    nz = [a for a in range(L.NA) if av[a]]
    nzc = [a for a in nz if a not in atom_out]
    nzg = [a for a in nz if a in atom_out]
    fail = L.failing_eqs(av)
    print(f'  iter {it}: nonzero atoms={len(nz)} (checks {len(nzc)}, gates {len(nzg)}) '
          f'failing={len(fail)} score={L.NEQ-len(fail)}')
    if it >= 1 and len(nz) == prev:
        break
    prev = len(nz)

print('\nnonzero CHECK atoms after forward eval:')
for a in nzc[:40]:
    print(f'   a{a:<6} neq={len(L.atom2eq.get(a,{})):<3} {L.atom_src[a][:120]}')
print(f'   ... total {len(nzc)}')
print('\nnonzero GATE atoms (cycles / non-exact divisions):')
for a in nzg[:20]:
    print(f'   a{a:<6} -> x_{atom_out[a][1]}  {L.atom_src[a][:100]}')
print(f'   ... total {len(nzg)}')
T.save(v, os.path.join(HERE, 'forward_state.json'))
