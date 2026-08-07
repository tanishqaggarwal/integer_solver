#!/usr/bin/env python3
"""Complete the triangulation by greedily designating free inputs."""
import pickle, os, collections, sys, heapq
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
freevars = list(pure)
used_atom = [False] * M.na
unk = [sum(1 for x in M.avars[i] if not produced[x]) for i in range(M.na)]
order = []
Q = collections.deque(i for i in range(M.na) if unk[i] == 1)


def drain():
    while Q:
        a = Q.popleft()
        if used_atom[a] or unk[a] != 1:
            continue
        miss = [x for x in M.avars[a] if not produced[x]]
        if len(miss) != 1 or miss[0] not in sol[a]:
            continue
        v = miss[0]
        produced[v] = True; used_atom[a] = True
        order.append((v, a))
        for b in var2atoms[v]:
            if not used_atom[b]:
                unk[b] -= 1
                if unk[b] == 1:
                    Q.append(b)


drain()
nnew = 0
while True:
    rem = [v for v in range(NV) if not produced[v]]
    if not rem:
        break
    # score: how many atoms with unk==2 contain v
    sc = collections.Counter()
    for a in range(M.na):
        if used_atom[a] or unk[a] != 2:
            continue
        for x in M.avars[a]:
            if not produced[x]:
                sc[x] += 1
    if not sc:
        v = rem[0]
    else:
        v = sc.most_common(1)[0][0]
    produced[v] = True
    freevars.append(v)
    nnew += 1
    for b in var2atoms[v]:
        if not used_atom[b]:
            unk[b] -= 1
            if unk[b] == 1:
                Q.append(b)
    drain()

print(f"extra free inputs designated: {nnew}")
print(f"total free inputs: {len(freevars)} (pure {len(pure)})")
print(f"DAG-produced: {len(order)}")
checks = [i for i in range(M.na) if not used_atom[i]]
print(f"check atoms: {len(checks)}")
pickle.dump({'order': order, 'checks': checks, 'freevars': freevars,
             'extra': nnew},
            open(os.path.join(HERE, f'orient3_{pol}.pkl'), 'wb'))
