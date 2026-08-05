#!/usr/bin/env python3
"""Downstream cone of the free inputs that feed the two violated checks."""
import pickle
from collections import defaultdict, deque
import heal_harness as H
import _om_parse as OP

D = pickle.load(open('_om_parsed2.pkl', 'rb')); astof = D['astof']
G = pickle.load(open('_om_dag.pkl', 'rb'))
usedatom = G['usedatom']; order = G['order']; checks = set(G['checks'])
avars = G['avars']; vat = G['vat']; free = set(G['free'])
pos = {v: i for i, v in enumerate(order)}

# children: var -> vars defined using it
kids = defaultdict(set)
for t, k in usedatom.items():
    for v in avars[k]:
        if v != t: kids[v].add(t)

def cone(seed):
    seen = set(seed); q = deque(seed)
    while q:
        v = q.popleft()
        for w in kids[v]:
            if w not in seen: seen.add(w); q.append(w)
    return seen

for seed in [[7068], [4432], [7068, 4432], [17325], [9413]]:
    C = cone(seed)
    ck = set()
    for v in C:
        for k in vat[v]:
            if k in checks: ck.add(k)
    print('cone(%s): %d vars, touches %d check atoms' % (seed, len(C), len(ck)))

# which free inputs feed x_2099 / x_19964 (the other side of the two checks)?
parents = defaultdict(set)
for t, k in usedatom.items():
    for v in avars[k]:
        if v != t: parents[t].add(v)
def anc(seed):
    seen = set(seed); q = deque(seed)
    while q:
        v = q.popleft()
        for w in parents[v]:
            if w not in seen: seen.add(w); q.append(w)
    return seen
for v in [2099, 19964, 7068, 4432]:
    A = anc([v])
    print('ancestors(x_%d): %d vars, %d of them free inputs' % (v, len(A), len(A & free)))
print()
print('x_17325 free?', 17325 in free, ' x_9413 free?', 9413 in free)
print('x_7068 free?', 7068 in free, ' x_4432 free?', 4432 in free)
print('x_28599 defined by:', usedatom.get(28599))
print('x_17499 defined by:', usedatom.get(17499))
# trace x_28599 back to a constant
v = 28599
for _ in range(12):
    k = usedatom.get(v)
    if k is None: print('  x_%d is FREE' % v); break
    print('  x_%d <- %s' % (v, k))
    ps = [w for w in avars[k] if w != v]
    if not ps: break
    v = ps[0]
