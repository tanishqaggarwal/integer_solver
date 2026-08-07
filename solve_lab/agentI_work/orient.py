#!/usr/bin/env python3
"""Structural triangulation: build a maximal acyclic definition DAG from the
pure inputs forward, leaving the rest as CHECK atoms."""
import pickle, os, collections, sys, json
from model import Model

HERE = os.path.dirname(os.path.abspath(__file__))
M = Model()
dag = pickle.load(open(os.path.join(HERE, 'dag.pkl'), 'rb'))
sol = dag['sol']          # per-atom: vars for which atom is unit-solvable
NV = 38748

var2atoms = collections.defaultdict(list)
for i, vs in enumerate(M.avars):
    for x in vs:
        var2atoms[x].append(i)

produced = [False] * NV        # variable value determined
definer = [None] * NV
# pure inputs
has_def = set()
for i, s in enumerate(sol):
    for v in s:
        has_def.add(v)
free = sorted(set(range(NV)) - has_def)
for v in free:
    produced[v] = True

# atom unknown-count
unk = [sum(1 for x in M.avars[i] if not produced[x]) for i in range(M.na)]
Q = collections.deque(i for i in range(M.na) if unk[i] == 1)
order = []
used_atom = [False] * M.na
while Q:
    a = Q.popleft()
    if used_atom[a] or unk[a] != 1:
        continue
    miss = [x for x in M.avars[a] if not produced[x]]
    if len(miss) != 1:
        continue
    v = miss[0]
    if v not in sol[a]:
        continue           # cannot solve exactly for it
    produced[v] = True
    definer[v] = a
    used_atom[a] = True
    order.append((v, a))
    for b in var2atoms[v]:
        if not used_atom[b]:
            unk[b] -= 1
            if unk[b] == 1:
                Q.append(b)

nprod = sum(produced)
print(f"pure inputs (no definer anywhere): {len(free)}")
print(f"variables produced by DAG: {len(order)}")
print(f"total produced: {nprod} / {NV}   unproduced: {NV - nprod}")
checks = [i for i in range(M.na) if not used_atom[i]]
print(f"check atoms (not used as definitions): {len(checks)}")
c = collections.Counter(sum(1 for x in M.avars[i] if not produced[x]) for i in checks)
print("checks by #unproduced vars:", sorted(c.items()))

pickle.dump({'free': free, 'order': order, 'definer': definer,
             'produced': produced, 'checks': checks},
            open(os.path.join(HERE, 'orient.pkl'), 'wb'))

# how many checks are trivially satisfied at the 39026 witness?
from model import load_assign
v = load_assign(os.path.join(HERE, '..', 'best', 'new_instance_partial_39026.json'))
bad = [i for i in checks if M.atom_val(i, v) != 0]
print("checks nonzero at witness:", len(bad), bad[:20])
