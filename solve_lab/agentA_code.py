#!/usr/bin/env python3
"""Find the CODE structure of the 78 message bits.
- forward cone (descendants) of the bits
- message-sensitive equations (checks the bits must satisfy)
- per-bit broken equation sets; analyze whether pairs/combos can cancel."""
import json, sys, re
from collections import defaultdict
from agentA_harness import (p, order, definer, gates, freeinp, backward_cone,
                            load_solution, forward, eqvars, eqcode, lines, NEQ, NVARS)

boolset = set(json.load(open('boolbits.json'))['boolvars'])
allS, freeS = backward_cone(35389)
allT, freeT = backward_cone(6671)
allbits = sorted((freeS | freeT) & boolset)
bitset = set(allbits)

# forward cone: gates whose vids include a bit or a downstream var
gdef_vids = {t: gates[definer[t]][2] for t in order}
consumers = defaultdict(list)   # var -> gate targets that use it
for t in order:
    for u in gdef_vids[t]:
        consumers[u].append(t)
fwd = set(bitset)
stack = list(bitset)
while stack:
    u = stack.pop()
    for t in consumers.get(u, ()):
        if t not in fwd:
            fwd.add(t); stack.append(t)
fwd_gates = fwd - bitset
print(f"forward cone of {len(bitset)} bits: {len(fwd_gates)} downstream gates")

# message-sensitive equations: eqs whose vars intersect fwd (bits or descendants)
msens = [i for i in range(NEQ) if eqvars[i] & fwd]
print(f"message-sensitive equations: {len(msens)}")

base = load_solution('best/new_instance_partial_39013.json')
forward(base)
F0 = set(i for i in msens if eval(eqcode[i], {'__builtins__': {}, 'v': base}) != 0)
print(f"of these, currently failing (core subset): {len(F0)} -> {sorted(F0)}")
sat = [i for i in msens if i not in F0]   # satisfied message-sensitive eqs = the checks
print(f"satisfied message-sensitive (the CHECK constraints): {len(sat)}")

# restricted forward eval over fwd cone only
fwd_order = [t for t in order if t in fwd_gates]
VAR = re.compile(r'x_(\d+)')
fwd_code = [compile(VAR.sub(r'v[\1]', gates[definer[t]][1]), '<r>', 'eval') for t in fwd_order]
ns = {'__builtins__': {}, 'v': base}
def pfwd(v):
    ns['v'] = v
    for c, t in zip(fwd_code, fwd_order):
        v[t] = eval(c, ns)

# per-bit break sets over the check constraints (sat)
satcode = [(i, eqcode[i]) for i in sat]
def broken(v):
    ns['v'] = v
    return set(i for i, c in satcode if eval(c, ns) != 0)

perbit = {}
for b in allbits:
    v = base[:]
    v[b] = 1 - base[b]
    pfwd(v)
    perbit[b] = broken(v)
allbroke = set().union(*perbit.values())
print(f"union of all single-bit-broken checks: {len(allbroke)}")
# these are the equations that constrain the message
sizes = sorted(len(perbit[b]) for b in allbits)
print(f"per-bit break-set sizes: min={sizes[0]} max={sizes[-1]} median={sizes[len(sizes)//2]}")

# overlap structure: build incidence bit x check
inc = {b: perbit[b] for b in allbits}
# how many checks are hit by exactly 1,2,... bits
hitcount = defaultdict(int)
for c in allbroke:
    n = sum(1 for b in allbits if c in inc[b])
    hitcount[n] += 1
print(f"checks hit by k bits: {dict(sorted(hitcount.items()))}")

json.dump({'allbits': allbits, 'msens': msens, 'sat_checks': sat, 'allbroke': sorted(allbroke),
           'perbit': {str(b): sorted(perbit[b]) for b in allbits}},
          open('agentA_code.json', 'w'))
print("saved agentA_code.json")
