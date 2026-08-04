#!/usr/bin/env python3
"""Probe the algebra: are checks GF(2)-parity-like? linear/quadratic in bits?
Test pair flips against symmetric-difference and second-difference hypotheses."""
import json, sys, re, random
from collections import defaultdict
from agentA_harness import (p, order, definer, gates, backward_cone,
                            load_solution, forward, eqcode, lines, eqvars, NEQ, NVARS)

boolset = set(json.load(open('boolbits.json'))['boolvars'])
cd = json.load(open('agentA_code.json'))
allbits = cd['allbits']; sat = cd['sat_checks']
perbit = {int(k): set(v) for k, v in cd['perbit'].items()}
fwd_gates_all = None

# rebuild fwd cone eval
allS, freeS = backward_cone(35389); allT, freeT = backward_cone(6671)
bitset = set(allbits)
gdef_vids = {t: gates[definer[t]][2] for t in order}
consumers = defaultdict(list)
for t in order:
    for u in gdef_vids[t]: consumers[u].append(t)
fwd = set(bitset); stack = list(bitset)
while stack:
    u = stack.pop()
    for t in consumers.get(u, ()):
        if t not in fwd: fwd.add(t); stack.append(t)
fwd_gates = fwd - bitset
fwd_order = [t for t in order if t in fwd_gates]
VAR = re.compile(r'x_(\d+)')
fwd_code = [compile(VAR.sub(r'v[\1]', gates[definer[t]][1]), '<r>', 'eval') for t in fwd_order]
base = load_solution('best/new_instance_partial_39013.json'); forward(base)
ns = {'__builtins__': {}, 'v': base}
def pfwd(v):
    ns['v'] = v
    for c, t in zip(fwd_code, fwd_order): v[t] = eval(c, ns)

satcode = [(i, eqcode[i]) for i in sat]
def broken(v):
    ns['v'] = v
    return set(i for i, c in satcode if eval(c, ns) != 0)

def flips(bits):
    v = base[:]
    for b in bits: v[b] = 1 - base[b]
    pfwd(v)
    return v

# --- Test GF(2) parity hypothesis on 25 random pairs ---
print("=== GF(2) parity test (pair-flip break == symdiff of single breaks?) ===")
random.seed(1)
match = 0; tot = 0
for _ in range(25):
    i, j = random.sample(allbits, 2)
    v = flips([i, j]); bij = broken(v)
    expect = perbit[i] ^ perbit[j]
    tot += 1
    if bij == expect: match += 1
    else:
        extra = bij - expect; missing = expect - bij
        if _ < 6: print(f"  ({i},{j}): actual {len(bij)} vs symdiff {len(expect)}; extra={len(extra)} missing={len(missing)}")
print(f"GF(2) symdiff matches: {match}/{tot}")

# --- Test S,T quadratic vs higher: 3rd order differences at baseline + random points ---
print("\n=== degree of S,T mod p in bits (3rd-order differences) ===")
def ST(v): return v[35389] % p, v[6671] % p
def d3(i, j, k, base_flips, which):
    def f(extra):
        v = flips(list(set(base_flips) ^ set(extra)))
        s, t = ST(v); return s if which == 'S' else t
    # 3rd finite difference
    return (f([i,j,k]) - f([i,j]) - f([i,k]) - f([j,k]) + f([i]) + f([j]) + f([k]) - f([])) % p
nz3S = 0; nz3T = 0; n3 = 0
sens = [b for b in allbits if (base[35389] != flips([b])[35389])]  # quick sensitive set
for _ in range(30):
    i, j, k = random.sample(allbits, 3)
    for bf in [[], random.sample(allbits, 3)]:
        n3 += 1
        if d3(i, j, k, bf, 'S') != 0: nz3S += 1
        if d3(i, j, k, bf, 'T') != 0: nz3T += 1
print(f"nonzero 3rd-order diffs: S={nz3S}/{n3}  T={nz3T}/{n3}  (0 => degree<=2)")

# --- Look at example check equations ---
print("\n=== example check equations ===")
# a check hit by 1 bit
by1 = [c for c in sat if sum(1 for b in allbits if c in perbit[b]) == 1]
bymany = [c for c in sat if sum(1 for b in allbits if c in perbit[b]) >= 41]
for label, lst in [('hit-by-1', by1[:3]), ('hit-by-many', bymany[:2])]:
    for c in lst:
        L = lines[c]
        print(f"  [{label}] eq {c}: len={len(L)} nvars={len(eqvars[c])}")
        print(f"     {L[:220]}")
