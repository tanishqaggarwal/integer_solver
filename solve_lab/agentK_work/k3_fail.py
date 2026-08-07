#!/usr/bin/env python3
"""K3: what are the 7 truly-failing equations made of? Which stage do they belong to?"""
import sys, os, json, collections
F = '/home/user/integer_solver/solve_lab/agentF_work'
sys.path.insert(0, F)
from fwd import Engine, NV, compile_node
from circ2 import vars_of

E = Engine()
FAIL = [12231, 12270, 12350, 14584, 18673, 22044, 29125]
d = json.load(open('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json'))
full = [0] * NV
for k, v in d.items():
    full[int(k[2:])] = int(v)

# where does each atom appear (which equations)
atomeq = collections.defaultdict(list)
for i, row in enumerate(E.eqrows):
    for k, a in row:
        atomeq[a].append(i)

for i in FAIL:
    row = E.eqrows[i]
    print('=== eq', i, 'natoms', len(row))
    for k, a in row:
        val = eval(compile(compile_node(E.atoms[a]), '<x>', 'eval'), {'v': full, '__builtins__': {}})
        if val:
            print('   coef', k, 'atom', a[:120], 'VAL!=0 vars', sorted(vars_of(E.atoms[a])))
    # summarize all atoms' variables
    vs = set()
    for k, a in row: vs |= vars_of(E.atoms[a])
    print('   total distinct vars', len(vs))

roles = json.load(open(F + '/stage_roles.json'))
# map: for each stage, which equations do its 3 checks live in?
print()
print('root 15298 wires -> equations containing atoms with those vars')
for st in ['15298', '30973', '24533']:
    r = roles[st][0]
    ws = set(r['out']) | set(r['inA']) | set(r['inB'])
    eqs = set()
    for a, node in E.atoms.items():
        if a in E.residx and vars_of(node) & ws:
            eqs |= set(atomeq[a])
    print(st, 'wires', sorted(ws), 'eqs', sorted(eqs)[:20], 'n', len(eqs))
