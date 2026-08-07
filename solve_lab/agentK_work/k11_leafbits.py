#!/usr/bin/env python3
"""K11: recover the liveness (selector) tree and its leaf bits from the atoms alone."""
import sys, os, json, collections, re
K = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, K)
from cascade import Cascade, NV, P
from parse import node_str

C = Cascade()
d = json.load(open('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json'))
full = [0] * NV
for k, val in d.items(): full[int(k[2:])] = int(val)
freeset = set(C.E.free)
vc = json.load(open(K + '/varclass.json'))
handles = set(vc['handles']); bools = set(vc['bools'])

# definition atoms of shape (xG-(xA*xB)) and (xG-(1-xB)) and (xG-xA)
prod = {}   # G -> (A,B)
alias = {}  # G -> A
notg = {}   # G -> A   meaning G = 1-A
for i, name in enumerate(C.atomnames):
    n = C.nodes[i]
    if n[0] != '-': continue
    l, r = n[1], n[2]
    if l[0] != 'v': continue
    G = l[1]
    if r[0] == '*' and r[1][0] == 'v' and r[2][0] == 'v':
        prod[G] = (r[1][1], r[2][1])
    elif r[0] == 'v':
        alias[G] = r[1]
    elif r[0] == '-' and r[1] == ('c', 1) and r[2][0] == 'v':
        notg[G] = r[2][1]
print('product defs', len(prod), 'alias defs', len(alias), 'not defs', len(notg))

# a var is "boolean" if it is a free bool, or an alias/not/product of boolean vars
isb = {}
def boolean(u, depth=0):
    if u in isb: return isb[u]
    if depth > 60: return False
    isb[u] = False
    if u in freeset:
        isb[u] = (u in bools); return isb[u]
    if u in alias: r = boolean(alias[u], depth + 1)
    elif u in notg: r = boolean(notg[u], depth + 1)
    elif u in prod:
        a, b = prod[u]; r = boolean(a, depth + 1) and boolean(b, depth + 1)
    else: r = False
    isb[u] = r
    return r

bset = [u for u in range(NV) if boolean(u)]
print('boolean-closure vars:', len(bset))
# leaves of the liveness tree: free boolean vars
leafbits = sorted(u for u in bset if u in freeset)
print('free boolean leaf-bit candidates:', len(leafbits))
on = [u for u in leafbits if full[u] == 1]
print('ON in deliverable:', len(on), on[:40])

# gate structure among boolean vars
gates = {u: prod[u] for u in bset if u in prod}
print('AND gates among boolean vars:', len(gates))
# resolve alias/not chains to base
def base(u, seen=None):
    while u in alias: u = alias[u]
    return u
# build support: leaf set under each boolean var
supp = {}
def support(u, depth=0):
    if u in supp: return supp[u]
    if depth > 60: return set()
    supp[u] = set()
    if u in freeset: s = {u}
    elif u in alias: s = support(alias[u], depth + 1)
    elif u in notg: s = support(notg[u], depth + 1)
    elif u in prod:
        a, b = prod[u]; s = support(a, depth + 1) | support(b, depth + 1)
    else: s = set()
    supp[u] = s
    return s
sizes = collections.Counter()
for u in bset: sizes[len(support(u))] += 1
print('support-size histogram (top):', sizes.most_common(20))
mx = max((len(support(u)), u) for u in bset)
print('largest boolean support:', mx[0], 'at x%d' % mx[1])
json.dump({'leafbits': leafbits, 'on_deliverable': on,
           'root_bool': mx[1], 'root_support': sorted(support(mx[1]))},
          open(K + '/leafbits.json', 'w'))
