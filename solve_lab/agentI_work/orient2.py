#!/usr/bin/env python3
"""Triangulate the residual circuit: seed = propagated values + pure inputs."""
import pickle, os, collections, sys
from model import Model
HERE = os.path.dirname(os.path.abspath(__file__))
M = Model()
dag = pickle.load(open(os.path.join(HERE, 'dag.pkl'), 'rb'))
sol = dag['sol']
NV = 38748
pol = sys.argv[1] if len(sys.argv) > 1 else 'wit'
val = pickle.load(open(os.path.join(HERE, f'prop_{pol}.pkl'), 'rb'))

var2atoms = collections.defaultdict(list)
for i, vs in enumerate(M.avars):
    for x in vs:
        var2atoms[x].append(i)

has_def = set()
for i, s in enumerate(sol):
    has_def |= set(s)
pure = sorted(set(range(NV)) - has_def)

produced = [False] * NV
for v in pure:
    produced[v] = True
for v in range(NV):
    if val[v] is not None:
        produced[v] = True
seed = sum(produced)
used_atom = [False] * M.na
unk = [sum(1 for x in M.avars[i] if not produced[x]) for i in range(M.na)]
Q = collections.deque(i for i in range(M.na) if unk[i] == 1)
order = []
while Q:
    a = Q.popleft()
    if used_atom[a] or unk[a] != 1:
        continue
    miss = [x for x in M.avars[a] if not produced[x]]
    if len(miss) != 1 or miss[0] not in sol[a]:
        continue
    v = miss[0]
    produced[v] = True; definer = a; used_atom[a] = True
    order.append((v, a))
    for b in var2atoms[v]:
        if not used_atom[b]:
            unk[b] -= 1
            if unk[b] == 1:
                Q.append(b)
print(f"seed produced {seed}; DAG produced {len(order)}; total {sum(produced)}/{NV}")
print(f"unproduced {NV - sum(produced)}")
checks = [i for i in range(M.na) if not used_atom[i]]
c = collections.Counter(sum(1 for x in M.avars[i] if not produced[x]) for i in checks)
print("check atoms:", len(checks), " by #unproduced:", sorted(c.items()))
pickle.dump({'order': order, 'produced': produced, 'checks': checks, 'pure': pure},
            open(os.path.join(HERE, f'orient2_{pol}.pkl'), 'wb'))
