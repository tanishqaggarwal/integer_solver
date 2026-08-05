#!/usr/bin/env python3
"""Custom augmenting-path re-matcher + exact forward re-eval.

Root cause of the 4 violated atoms: greedy propagation used *combination* atoms
(linear combos of gate residuals, nv>4 — redundant checks) to DEFINE variables,
in the wrong topological order. That pinned x_9770 (via 741-mon atom 40782) and
x_3183 before their true primitive sum-gates (27973, 27978) could fire, leaving
those gates violated -> which makes 41470 (contains 27973) and 45004 (contains
27978) violated too.

Fix: re-orient the variable<->gate matching so each variable is defined by a
PRIMITIVE gate (nv<=4). Combination atoms are redundant (implied by primitives)
and may be left unmatched. We find augmenting paths for the violated primitive
gates, then forward-evaluate the affected cone exactly over Z."""
import json, sys, time
from collections import defaultdict, deque

t0 = time.time()
NV_PRIM = 4  # atom is a primitive gate iff it touches <= this many variables

# ---- load atoms as list[(monomial_tuple, coef)] ----
atoms = []
avars = []
for line in open('atoms/poly_atoms.jsonl'):
    d = json.loads(line)
    poly = [(tuple(m), c) for m, c in d['poly']]
    atoms.append(poly)
    vs = set()
    for m, c in poly: vs.update(m)
    avars.append(vs)
NA = len(atoms)
NVARS = 1 + max(max(vs) if vs else 0 for vs in avars)
print(f"loaded {NA} atoms, NVARS={NVARS} ({time.time()-t0:.0f}s)", flush=True)

# ---- current matching from recorded provenance ----
eo = json.load(open('eval_order.json')); prov = eo['prov']
owner = [None]*NVARS          # var -> defining atom
defines = [None]*NA           # atom -> var it defines
for v, p in enumerate(prov):
    if p and p[0] >= 0:
        owner[v] = p[0]; defines[p[0]] = v

is_comb = [len(avars[a]) > NV_PRIM for a in range(NA)]

# candidate output vars of an atom: vars that appear ONLY in degree-1 monomials
def out_cands(a):
    bad = set(); lin = set()
    for m, c in atoms[a]:
        if len(m) == 1: lin.add(m[0])
        else:
            for x in m: bad.add(x)
    return [v for v in lin if v not in bad]
outc = [out_cands(a) for a in range(NA)]

# ---- find violated primitive atoms in current best solution ----
best = json.load(open('best/best_partial_39013.json'))
val = [0]*NVARS
for k, x in best.items():
    val[int(k[2:])] = x
def atom_val(a):
    s = 0
    for m, c in atoms[a]:
        t = c
        for x in m: t *= val[x]
        s += t
    return s
violated_prim = [a for a in range(NA) if not is_comb[a] and defines[a] is None and atom_val(a) != 0]
# also primitive atoms that ARE matched but still evaluate nonzero (shouldn't be)
violated_matched = [a for a in range(NA) if not is_comb[a] and defines[a] is not None and atom_val(a) != 0]
print(f"violated primitive atoms (unmatched): {violated_prim}", flush=True)
print(f"violated primitive atoms (matched!): {violated_matched}", flush=True)

# ---- augmenting-path re-match (each var defined by a primitive gate) ----
def augment(start_atom):
    """Find alternating path so start_atom claims a var; primitive atoms it
    displaces must be re-matched; combination atoms / free vars terminate."""
    stackpath = []                 # list of (atom, var) edges chosen
    visited_atoms = set([start_atom])
    def dfs(a):
        for v in outc[a]:
            o = owner[v]
            if o == a: continue
            if o is None:          # free var -> terminate
                stackpath.append((a, v)); return True
            if is_comb[o]:         # displace a redundant combination atom
                stackpath.append((a, v)); defines[o] = None; return True
            if o in visited_atoms: continue
            visited_atoms.add(o)
            stackpath.append((a, v))
            if dfs(o): return True
            stackpath.pop()
        return False
    if dfs(start_atom):
        for a, v in stackpath:
            owner[v] = a; defines[a] = v
        return True
    return False

changed_vars = set()
for a in violated_prim:
    old = list(owner)
    if augment(a):
        for v in range(NVARS):
            if owner[v] != old[v]: changed_vars.add(v)
        print(f"  augmented atom {a} -> now defines x_{defines[a]}", flush=True)
    else:
        print(f"  !! could not augment atom {a}", flush=True)
print(f"changed variable definitions: {sorted(changed_vars)} ({time.time()-t0:.0f}s)", flush=True)

# ---- forward re-evaluation of the affected cone, exact over Z ----
# dependency: var v (defined by atom a) depends on the other vars of a.
users = defaultdict(list)   # var -> list of vars whose defining atom uses it
for v in range(NVARS):
    a = owner[v]
    if a is None: continue
    for x in avars[a]:
        if x != v: users[x].append(v)

# forward closure of changed vars
dirty = set(changed_vars)
dq = deque(changed_vars)
while dq:
    x = dq.popleft()
    for u in users[x]:
        if u not in dirty:
            dirty.add(u); dq.append(u)
print(f"cone to recompute: {len(dirty)} vars", flush=True)

# topological order within dirty (Kahn over dependency edges among dirty vars)
indeg = {v: 0 for v in dirty}
adj = defaultdict(list)
for v in dirty:
    a = owner[v]
    for x in avars[a]:
        if x != v and x in dirty:
            adj[x].append(v); indeg[v] += 1
q = deque([v for v in dirty if indeg[v] == 0])
topo = []
while q:
    v = q.popleft(); topo.append(v)
    for u in adj[v]:
        indeg[u] -= 1
        if indeg[u] == 0: q.append(u)
if len(topo) != len(dirty):
    print(f"  !! cycle in cone: topo {len(topo)} < dirty {len(dirty)} — cyclic vars remain", flush=True)

def solve_for(a, v):
    """Return exact integer value of v s.t. atom a == 0, or None if not divisible."""
    coef_v = 0; rest = 0
    for m, c in atoms[a]:
        if m == (v,):
            coef_v += c
        else:
            t = c
            for x in m: t *= val[x]
            rest += t
    if coef_v == 0: return None
    if (-rest) % coef_v != 0: return None
    return (-rest) // coef_v

fail = 0
for v in topo:
    nvv = solve_for(owner[v], v)
    if nvv is None:
        fail += 1
        if fail <= 10: print(f"  non-divisible solve for x_{v} via atom {owner[v]}", flush=True)
        continue
    val[v] = nvv
print(f"recomputed {len(topo)} vars, {fail} non-divisible ({time.time()-t0:.0f}s)", flush=True)

# ---- full violation check over Z ----
viol = []
for a in range(NA):
    if atom_val(a) != 0: viol.append(a)
print(f"TOTAL violated atoms after re-solve: {len(viol)}", flush=True)
print("  first 30:", viol[:30], flush=True)
json.dump({f"x_{i}": val[i] for i in range(NVARS)}, open('cand_rematch.json', 'w'))
print("wrote cand_rematch.json", flush=True)
