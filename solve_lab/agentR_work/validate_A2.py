#!/usr/bin/env python3
"""Is the (2)/(3) absence real, or just because the deliverable is a PARTIAL assignment?
pins.json names, for each pin variable, the two WIRES that should carry its point's coordinates.
For the two pins that ARE set to 1 in the deliverable, check those exact wires."""
import sys, json
sys.path.insert(0, '.')
import model
from model import P, S
D = json.load(open('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json'))
dv = {}
for k, v in D.items():
    dv[int(k[2:]) if k.startswith('x_') else int(k)] = int(v)
pins = json.load(open('/home/user/integer_solver/solve_lab/agentF_work/pins.json'))
lad = [(int(a), int(b)) for _, a, b in json.load(open('ladder.json'))['ladder']]
LX = {x: i for i, (x, y) in enumerate(lad)}
print('deliverable specifies %d of 38748 variables' % len(dv))
for v in (24601, 2081):
    e = pins[str(v)]
    print('\npin x%d = %s' % (v, dv.get(v, 0)))
    for slot, (w, val) in enumerate(e):
        w = int(w); val = int(val)
        present = w in dv
        li = LX.get((val + S) % P)
        print('   its wire x%-6d (ladder idx %s):  set in deliverable? %s   value there = %s'
              % (w, li, present, (str(dv[w])[:38] + '...') if present else 'UNSET(=0)'))
        print('        pins.json says that wire should hold %s...' % str(val)[:38])
        if present: print('        MATCH' if dv[w] == val else '        MISMATCH')
# how many of the 512 coordinate wires named by pins.json are set at all?
wires = []
for v, e in pins.items():
    for w, val in e:
        wires.append((int(w), int(val)))
setn = sum(1 for w, _ in wires if w in dv)
match = sum(1 for w, val in wires if dv.get(w) == val)
print('\ncoordinate wires named by pins.json: %d ; set in deliverable: %d ; holding the named value: %d'
      % (len(wires), setn, match))
json.dump({'n_specified': len(dv), 'coord_wires': len(wires), 'set': setn, 'match': match},
          open('runs/validate_A2.json', 'w'), indent=1)
