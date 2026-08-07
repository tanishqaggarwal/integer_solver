#!/usr/bin/env python3
"""K25: role-based variable classification straight from atom shapes.
  BOOL   : free var u with an idempotency atom (xu*(xu-1))
  HANDLE : free var whose every atom coefficient is divisible by p
  WIRE   : every other free var (leaf wires, slot wires, stage outputs)"""
import sys, os, json, re, collections
K = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, K)
F = '/home/user/integer_solver/solve_lab/agentF_work'
sys.path.insert(0, F)
from cascadep import CascadeP, NV, P

C = CascadeP()
freeset = set(C.E.free)
boolv = set()
for nm in C.names:
    m = re.fullmatch(r'\(x(\d+)\*\(x(\d+)-1\)\)', nm)
    if m and m.group(1) == m.group(2): boolv.add(int(m.group(1)))
print('idempotency atoms -> boolean vars:', len(boolv), 'of which free:', len(boolv & freeset))

d = json.load(open('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json'))
full = [0] * NV
for k, val in d.items(): full[int(k[2:])] = int(val)

# handle detection: linear coefficient of u in every atom containing u is == 0 mod p
handles = []
for u in sorted(freeset):
    cs = []
    for i in C.var2atoms[u]:
        old = full[u]
        full[u] = old + 1; c1 = C.ev(i, full)
        full[u] = old; c0 = C.ev(i, full)
        cs.append((c1 - c0) % P)
    if cs and all(c == 0 for c in cs): handles.append(u)
handles = [u for u in handles if u not in boolv]
print('handles (free, coeff==0 mod p everywhere):', len(handles))
bools = sorted(boolv & freeset)
wires = sorted(freeset - set(handles) - set(bools))
print('free split: bools %d, handles %d, wires %d, total %d (free=%d)'
      % (len(bools), len(handles), len(wires), len(bools) + len(handles) + len(wires), len(freeset)))
D = json.load(open(K + '/points.json'))
leafsel = set(l['sel'] for l in D['leaves'])
print('leaf selectors that are free bools:', len(leafsel & set(bools)), 'of 256')
print('non-leaf free bools:', len(set(bools) - leafsel))
json.dump({'bools': bools, 'handles': handles, 'wires': wires,
           'leafsel': sorted(leafsel), 'otherbools': sorted(set(bools) - leafsel)},
          open(K + '/varclass2.json', 'w'))
