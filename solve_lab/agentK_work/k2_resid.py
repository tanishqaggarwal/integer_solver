#!/usr/bin/env python3
"""K2: with the deliverable's FULL assignment, which residual atoms are nonzero and what do they demand?"""
import sys, os, json, pickle, time, collections
F = '/home/user/integer_solver/solve_lab/agentF_work'
sys.path.insert(0, F)
from fwd import Engine, NV
from parse import node_str
from circ2 import vars_of

p = 115792089237316195423570985008687907853269984665640564039457584007908834671663
K = 97553848499418123410591666447050222001188385549510401465815187079080512838891

E = Engine()
d = json.load(open('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json'))
full = [0] * NV
for k, v in d.items():
    full[int(k[2:])] = int(v)

# evaluate residual atoms directly at the deliverable assignment (no forward overwrite)
r = [0] * len(E.res)
exec(E.rprog, {'v': full, 'r': r, '__builtins__': {}})
nz = [i for i, x in enumerate(r) if x]
print('nonzero residual atoms at deliverable:', len(nz))
for i in nz:
    print('   ', E.res[i], ' value mod p ==0?', r[i] % p == 0)
bad = E.score(r)
print('failing equations:', len(bad), bad)

# also verify definition atoms are zero
defbad = []
for a in E.order:
    c = E.cls[a]
    val = eval(compile(__import__('fwd').compile_node(E.atoms[a]), '<x>', 'eval'), {'v': full, '__builtins__': {}})
    if val: defbad.append(a)
print('nonzero definition atoms:', len(defbad), defbad[:10])

roles = json.load(open(F + '/stage_roles.json'))
r0 = roles['15298'][0]
print('root wires', r0)
for nm in ('out', 'inA', 'inB'):
    print(nm, [full[w] % p for w in r0[nm]])
