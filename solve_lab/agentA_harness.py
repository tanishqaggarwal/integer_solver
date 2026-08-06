#!/usr/bin/env python3
"""Reusable harness for the boolean-message algebraic attack.
Builds the gate DAG topo order (like forward_construct.py), provides forward-eval,
backward free-input cone computation, and loads the best baseline solution.
"""
import json, re, ast, sys
from collections import defaultdict, deque
from propagate import load_atoms, atom_vars, NVARS
sys.setrecursionlimit(1000000)

p = 2**256 - 2**32 - 977
hc = json.load(open('huge_consts.json')); C1 = int(hc['C1']); C2 = int(hc['C2'])

A = load_atoms()
gates = []
with open('atoms/gates.jsonl') as f:
    for line in f:
        d = json.loads(line); gates.append((d['t'], d['rhs'], tuple(d['vids'])))

# ---- single-var atom pins (same as forward_construct) ----
pin_from_atom = {}
for pp in A:
    vs = atom_vars(pp)
    if len(vs) == 1:
        v = next(iter(vs)); c0 = pp.get((), 0); c1 = pp.get((v,), 0); c2 = pp.get((v, v), 0)
        if c2 == 0 and c1 != 0 and (-c0) % c1 == 0:
            pin_from_atom.setdefault(v, (-c0) // c1)

gate_out = set(t for t, _, _ in gates)
freeinp = set(v for v in range(NVARS) if v not in gate_out)
override = {24601: 1, 2081: 1, 30213: C2, 22162: C1, 24468: C1, 18956: C2}

# ---- build topo order over gates ----
cand = defaultdict(list)
for gi, (t, rhs, vids) in enumerate(gates): cand[t].append(gi)
targets = set(cand)

def build_order():
    val = [0]*NVARS; pinned = [False]*NVARS
    for v, x in pin_from_atom.items():
        if not pinned[v]: val[v] = x; pinned[v] = True
    for v, x in override.items(): val[v] = x; pinned[v] = True
    ready = [False]*NVARS
    for v in range(NVARS):
        if v not in targets or v in freeinp or pinned[v]: ready[v] = True
    gu = [0]*len(gates); using = defaultdict(list)
    for gi, (t, rhs, vids) in enumerate(gates):
        u = 0
        for v in vids:
            if not ready[v]: u += 1
            using[v].append(gi)
        gu[gi] = u
    definer = {}; order = []
    q = deque(gi for gi in range(len(gates)) if gu[gi] == 0)
    while q:
        gi = q.popleft(); t, rhs, vids = gates[gi]
        if ready[t]: continue
        definer[t] = gi; order.append(t); ready[t] = True
        for gj in using[t]:
            gu[gj] -= 1
            if gu[gj] == 0: q.append(gj)
    return val, pinned, definer, order

val0, pinned, definer, order = build_order()
VAR = re.compile(r'x_(\d+)')
gcode = [compile(VAR.sub(r'v[\1]', gates[definer[order[k]]][1]), '<r>', 'eval') for k in range(len(order))]

# free-input ancestors of each gate target
anc = defaultdict(set)
for v in freeinp: anc[v] = {v}
for k, t in enumerate(order):
    _, rhs, vids = gates[definer[t]]
    s = set()
    for u in vids: s |= anc[u]
    anc[t] = s

# ---- equations ----
lines = [L for L in open('../EQUATIONS.txt').read().split('\n') if L.strip()]
eqcode = [compile(VAR.sub(r'v[\1]', L.rsplit('=', 1)[0]), '<e>', 'eval') for L in lines]
eqvars = [set(int(m) for m in VAR.findall(L)) for L in lines]
NEQ = len(lines)
ns = {'__builtins__': {}}

def forward(val):
    ns['v'] = val
    for k, t in enumerate(order): val[t] = eval(gcode[k], ns)
    return val

def eval_fails(val, subset=None):
    ns['v'] = val
    idxs = range(NEQ) if subset is None else subset
    return [i for i in idxs if eval(eqcode[i], ns) != 0]

def load_solution(path):
    d = json.load(open(path))
    val = [0]*NVARS
    for k, x in d.items():
        idx = int(k[2:]) if k.startswith('x_') else int(k)
        val[idx] = int(x)
    return val

# backward cone: all gate/free vars that S depends on (transitively)
gdef_vids = {}
for k, t in enumerate(order):
    gdef_vids[t] = gates[definer[t]][2]

def backward_cone(root):
    """Return (all_vars, free_inputs) in the backward cone of root var."""
    seen = set(); stack = [root]
    while stack:
        u = stack.pop()
        if u in seen: continue
        seen.add(u)
        for w in gdef_vids.get(u, ()):
            if w not in seen: stack.append(w)
    frees = set(u for u in seen if u in freeinp)
    return seen, frees

if __name__ == '__main__':
    boolset = set(json.load(open('boolbits.json'))['boolvars'])
    print(f"NVARS={NVARS} gates={len(gates)} order={len(order)} freeinp={len(freeinp)} boolvars={len(boolset)}")
    val = load_solution('best/new_instance_partial_39013.json')
    # check forward-eval consistency
    val2 = val[:]
    forward(val2)
    diff = [t for t in order if val2[t] != val[t]]
    print(f"forward-eval changed {len(diff)} gate targets vs stored solution")
    F = eval_fails(val2)
    print(f"after forward-eval: {NEQ-len(F)}/{NEQ} satisfied ({len(F)} fail)")
    print(f"S=x_35389 mod p = {val2[35389]%p}")
    print(f"T=x_6671  mod p = {val2[6671]%p}")
    for root, name in [(35389, 'S'), (6671, 'T')]:
        allv, frees = backward_cone(root)
        bfree = frees & boolset
        sfree = frees - boolset
        print(f"{name}=x_{root}: cone {len(allv)} vars, {len(frees)} free inputs; {len(bfree)} boolean, {len(sfree)} slack")
