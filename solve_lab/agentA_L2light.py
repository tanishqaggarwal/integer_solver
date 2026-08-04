#!/usr/bin/env python3
"""Minimal-perturbation L2≡0 mod 6672769 config (does NOT force S,T=0). Shift x_3558 to
a root mod M via a tiny x_16742 nudge δ<M, align x_19083 (x_16787+=δ), cascade-heal.
Report total + non-core wiring (the 20 core fail in wire=p since S,T!=0; wire=1 fixes them)."""
import json, re, sys
from collections import deque
from agentA_harness import (p, order, definer, gates, freeinp, backward_cone,
                            load_solution, forward, eqcode, eqvars, lines, NEQ, NVARS)
sys.setrecursionlimit(1000000)
M = 6672769
base = load_solution('best/new_instance_partial_39013.json'); forward(base)
gvids = {t: gates[definer[t]][2] for t in order}
CORE = {2071, 4573, 7123, 7469, 11854, 13660, 15299, 16622, 17726, 21382, 22093, 25480,
        25539, 28653, 29437, 31061, 32894, 32916, 34517, 34892}
def count(v):
    ns = {'__builtins__': {}, 'v': v}
    return set(i for i in range(NEQ) if eval(eqcode[i], ns) != 0)

# roots for x_3558 mod M
roots = [2783706, 5883594]
y0 = base[3558]  # baseline x_3558
deltas = [((y0 - r) % M, r) for r in roots]
deltas.sort()
dmin, root = deltas[0]
# choose sign to minimize |shift|: x_16742 += d makes x_3558 -= d ; want x_3558 ≡ root
# x_3558_new = y0 - d ≡ root  => d ≡ y0-root mod M ; pick smallest magnitude in (-M/2,M/2]
d = (y0 - root) % M
if d > M//2: d -= M
print(f"baseline x_3558 mod M = {y0%M}; target root = {root}; x_16742 nudge δ = {d} (|δ|<M)")

crit = {16742, 12186, 24908}
for r in [35389, 6671]:
    _, fr = backward_cone(r); crit |= set(fr)
QUOT = {30317, 2936, 5146}
protect = set(crit) | QUOT

from collections import Counter
allvarcount = Counter()
for i in range(NEQ):
    for x in eqvars[i]: allvarcount[x] += 1
const_pinned = set()
_bp = re.compile(r'\(x_(\d+)\)\s*-\s*\(?-?\d{20,}')
_bp2 = re.compile(r'-?\d{20,}\)?\s*-\s*\(x_(\d+)\)')
for i in range(NEQ):
    L = lines[i]
    if len(L) < 4000:
        for m in _bp.finditer(L): const_pinned.add(int(m.group(1)))
        for m in _bp2.finditer(L): const_pinned.add(int(m.group(1)))

_knob = {}
def free_cone(root_):
    seen=set(); st=[root_]
    while st:
        u=st.pop()
        if u in seen: continue
        seen.add(u)
        for w in gvids.get(u, ()):
            if w not in seen: st.append(w)
    return [u for u in seen if u in freeinp]
def additive_knob(v, gate):
    if gate in _knob: return _knob[gate]
    g0 = v[gate]; cands = []
    for w in free_cone(gate):
        if w in protect or w in const_pinned: continue
        old = v[w]; v[w] = old+1; forward(v); dd = v[gate]-g0; v[w] = old
        if dd in (1, -1): cands.append((allvarcount[w], w, dd))
    forward(v); cands.sort()
    res = (cands[0][1], cands[0][2]) if cands else (None, None)
    _knob[gate] = res; return res
def partners_in(eq_i, f):
    L = lines[eq_i]; out = set()
    for m in re.finditer(r'\(x_' + str(f) + r'\)\s*-\s*\(x_(\d+)\)', L): out.add(int(m.group(1)))
    for m in re.finditer(r'\(x_(\d+)\)\s*-\s*\(x_' + str(f) + r'\)', L): out.add(int(m.group(1)))
    return out
_v2e = {}
def eqs_with(var):
    if var not in _v2e: _v2e[var] = [i for i in range(NEQ) if var in eqvars[i]]
    return _v2e[var]

# apply minimal nudge: x_16742 += d ; keep x_16742 = x_19083 by aligning x_19083 (x_16787 += d)
v = base[:]
v[16742] = base[16742] + d
forward(v)
# align x_19083 = x_16742
def set_quot(v):
    pass  # in wire=p, quotient handles need S,T=0 which we don't have; leave core to wire=1
cur = count(v)
print(f"after x_16742 nudge (pre-align): {NEQ-len(cur)}/{NEQ} ({len(cur)} fail)")

assigned = {16742: v[16742]}
queue = deque([16742])
aligned_gates = set(); best = (len(cur), [x for x in v]); rounds = 0
while queue and rounds < 20000:
    rounds += 1
    f = queue.popleft(); fval = v[f]
    partners = set()
    for i in eqs_with(f): partners |= partners_in(i, f)
    for G in partners:
        if G in freeinp or G in aligned_gates: continue
        if v[G] == fval: continue
        w, dd = additive_knob(v, G)
        if w is None: continue
        old = v[w]; v[w] = old + dd * (fval - v[G]); forward(v)
        if v[G] != fval:
            v[w] = old; forward(v); continue
        aligned_gates.add(G); cur = count(v)
        if w not in assigned: assigned[w] = v[w]; queue.append(w)
        if len(cur) < best[0]: best = (len(cur), [x for x in v])
    if not cur: break
v = best[1]; cur = count(v)
noncore = cur - CORE
print(f"\nFINAL light config: {NEQ-len(cur)}/{NEQ} total fail; non-core fail={len(noncore)}: {sorted(noncore)[:25]}")
print(f"L2 = x_25739 mod {M} = {v[25739]%M}  (want 0)")
print(f"x_3558 mod M = {v[3558]%M} (want {root})  S,T mod p zero? {v[35389]%p==0},{v[6671]%p==0}")
print(f"assigned {len(assigned)} free inputs; core still failing (expected in wire=p): {len(cur & CORE)}")
if v[25739] % M == 0:
    json.dump({f"x_{i}": v[i] for i in range(NVARS) if v[i] != 0}, open('agentA_L2zero_light.json', 'w'))
    print("SAVED agentA_L2zero_light.json (L2≡0 mod 6672769, S,T!=0, minimal perturbation)")
