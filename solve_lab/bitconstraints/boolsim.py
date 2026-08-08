#!/usr/bin/env python3
"""Pass 8: pure boolean-layer simulator.

Extract the sub-circuit built only from
    copy  x - y,  not  1 - x - y,  and  x*y - z,  pin1 1 - x,  pin0 x
restricted to wires reachable from the 256 selector bits, evaluate it for many
bit vectors, and check every pin that lands on a reachable wire.  A pin that is
violated for some bit vectors and satisfied for others is a genuine pure-bit
constraint; we then characterise it exactly.
"""
import os, json, random, itertools
from collections import Counter, defaultdict
from scan import load, support, degree, fmt

HERE = os.path.dirname(os.path.abspath(__file__))
D, chain, BITS = load()
atoms = D['atoms']
atom_eqs = defaultdict(list)
for i, terms in enumerate(D['eq_terms']):
    for c, aid in terms:
        atom_eqs[aid].append(i)

copies, nots, ands, pin1, pin0 = [], [], [], {}, {}
for i, a in enumerate(atoms):
    d = dict(a); s = support(a); dg = degree(a)
    if len(a) == 1 and len(s) == 1 and dg == 1:
        pin0.setdefault(next(iter(s)), i)
    elif len(a) == 2 and len(s) == 1 and dg == 1 and () in d:
        v = next(iter(s))
        if d[()] == 1 and d[(v,)] == -1:
            pin1.setdefault(v, i)
    elif len(a) == 2 and len(s) == 2 and dg == 1:
        u, v = sorted(s)
        if d[(u,)] == 1 and d[(v,)] == -1:
            copies.append((i, u, v))
    elif len(a) == 3 and len(s) == 2 and dg == 1 and () in d:
        u, v = sorted(s)
        if d[()] == 1 and d[(u,)] == -1 and d[(v,)] == -1:
            nots.append((i, u, v))
    elif len(a) == 2 and len(s) == 3 and dg == 2:
        q = [m for m in d if len(m) == 2]; l = [m for m in d if len(m) == 1]
        if len(q) == 1 and len(l) == 1 and d[q[0]] == -d[l[0]] and q[0][0] != q[0][1]:
            ands.append((i, q[0][0], q[0][1], l[0][0]))

# ---- schedule ----
adj = defaultdict(list)          # input wire -> list of (kind, args..., out, atom)
for i, u, v in copies:
    adj[u].append(('copy', u, v, i)); adj[v].append(('copy', v, u, i))
for i, u, v in nots:
    adj[u].append(('not', u, v, i)); adj[v].append(('not', v, u, i))
for i, u, v, w in ands:
    adj[u].append(('and', u, v, w, i)); adj[v].append(('and', u, v, w, i))

sched = []           # ordered list of ops
known = set(chain)
frontier = list(chain)
prov = {b: [] for b in chain}
while frontier:
    nxt = []
    for x in frontier:
        for g in adj[x]:
            if g[0] in ('copy', 'not'):
                _, a, b, i = g
                if b not in known:
                    known.add(b); sched.append(g); prov[b] = prov[a] + [i]; nxt.append(b)
            else:
                _, u, v, w, i = g
                if w not in known and u in known and v in known:
                    known.add(w); sched.append(g); prov[w] = prov[u] + prov[v] + [i]
                    nxt.append(w)
    frontier = nxt
print(f"boolean wires reachable from the 256 selectors: {len(known)} "
      f"(ops: {Counter(g[0] for g in sched)})")


def simulate(bv):
    val = {}
    for b, x in zip(chain, bv):
        val[b] = x
    for g in sched:
        if g[0] == 'copy':
            val[g[2]] = val[g[1]]
        elif g[0] == 'not':
            val[g[2]] = 1 - val[g[1]]
        else:
            val[g[3]] = val[g[1]] * val[g[2]]
    return val


checks = [(w, 1, pin1[w]) for w in pin1 if w in known] + \
         [(w, 0, pin0[w]) for w in pin0 if w in known]
print(f"pins landing on reachable boolean wires: {len(checks)} "
      f"({sum(1 for c in checks if c[1]==1)} to 1, "
      f"{sum(1 for c in checks if c[1]==0)} to 0)")

rng = random.Random(7)
tests = [('zeros', [0]*256), ('ones', [1]*256)]
for k in range(20):
    tests.append((f'rand{k}', [rng.randrange(2) for _ in range(256)]))
for j in (0, 1, 7, 42, 128, 255):
    tests.append((f'e{j}', [1 if t == j else 0 for t in range(256)]))

viol = defaultdict(list)
for label, bv in tests:
    val = simulate(bv)
    for w, target, aid in checks:
        if val[w] != target:
            viol[(w, target, aid)].append(label)

print(f"\npins violated by at least one test vector: {len(viol)}")
records = []
for (w, target, aid), labels in sorted(viol.items()):
    print(f"  x_{w} pinned to {target} (atom#{aid}, eqs {sorted(set(atom_eqs[aid]))[:8]})"
          f" violated by {labels[:8]}{'...' if len(labels)>8 else ''} "
          f"({len(labels)}/{len(tests)})")

# ---- characterise: for each violated pin, find its exact boolean function ----
def fn_support(w):
    """which selectors actually influence w (single-bit flip test)"""
    base = [0]*256
    v0 = simulate(base)[w]
    dep = []
    for j in range(256):
        bv = list(base); bv[j] = 1
        if simulate(bv)[w] != v0:
            dep.append(j)
    return v0, dep


summary = []
for (w, target, aid), labels in sorted(viol.items()):
    v0, dep = fn_support(w)
    # test OR hypothesis: w == OR(dep)?
    ok = True
    for _ in range(60):
        bv = [0]*256
        for j in rng.sample(dep, rng.randrange(0, min(len(dep), 8)+1)):
            bv[j] = 1
        for j in rng.sample(range(256), 20):
            bv[j] = rng.randrange(2)
        got = simulate(bv)[w]
        want = 1 if any(bv[j] for j in dep) else 0
        if got != want:
            ok = False; break
    print(f"  -> x_{w}: value at all-zeros = {v0}; sensitive to {len(dep)} bits; "
          f"OR-hypothesis {'HOLDS' if ok else 'FAILS'}")
    summary.append({'wire': w, 'pin_value': target, 'pin_atom': aid,
                    'pin_equations': sorted(set(atom_eqs[aid])),
                    'value_at_all_zero': v0,
                    'sensitive_bits': dep,
                    'is_OR_of_sensitive_bits': ok,
                    'gate_atoms': prov[w]})

json.dump({'n_bool_wires': len(known), 'n_pins_checked': len(checks),
           'violated': summary},
          open(os.path.join(HERE, 'boolsim.json'), 'w'))
print("wrote boolsim.json")
