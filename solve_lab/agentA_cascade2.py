#!/usr/bin/env python3
"""Precise partner-alignment cascade. Track dirty free inputs; for each, find its gate
partner G in broken difference-checks (regex (x_f)-(x_G) / (x_G)-(x_f)); align G's value
to v[f] via G's live additive free knob; enqueue that knob. Recurse to fixpoint."""
import json, re, sys
from collections import defaultdict, deque
from agentA_harness import (p, order, definer, gates, freeinp, backward_cone,
                            load_solution, forward, eqcode, eqvars, lines, NEQ, NVARS)
sys.setrecursionlimit(1000000)
base = load_solution('best/new_instance_partial_39013.json'); forward(base)
gvids = {t: gates[definer[t]][2] for t in order}

crit = {16742, 14853, 12186, 24908}
for r in [35389, 6671]:
    _, fr = backward_cone(r); crit |= set(fr)
QUOT = {30317, 2936, 5146}
protect = set(crit) | QUOT

def set_quot(v):
    if v[11150] % p == 0: v[30317] = -(v[11150])//p
    if (537773*v[37758]) % p == 0: v[2936] = (537773*v[37758])//p
    if v[25739] % (6672769*p) == 0: v[5146] = v[25739]//(6672769*p)
def count(v):
    ns = {'__builtins__': {}, 'v': v}
    return set(i for i in range(NEQ) if eval(eqcode[i], ns) != 0)

# additive-knob cache
_knob = {}
def free_cone(root, limit=600):
    seen=set(); st=[root]
    while st and len(seen)<limit:
        u=st.pop()
        if u in seen: continue
        seen.add(u)
        for w in gvids.get(u, ()):
            if w not in seen: st.append(w)
    return [u for u in seen if u in freeinp]
def additive_knob(v, gate):
    if gate in _knob: return _knob[gate]
    g0 = v[gate]; res = (None, None)
    for w in free_cone(gate):
        if w in protect: continue
        old = v[w]; v[w] = old+1; forward(v); d = v[gate]-g0; v[w] = old
        if d in (1, -1): res = (w, d); break
    forward(v); _knob[gate] = res; return res

# var -> equations (only scan lazily)
_var2eq = {}
def eqs_with(var):
    if var not in _var2eq:
        _var2eq[var] = [i for i in range(NEQ) if var in eqvars[i]]
    return _var2eq[var]

def partners_in(eq_i, f):
    """gate partners G such that (x_f)-(x_G) or (x_G)-(x_f) appears in eq_i."""
    L = lines[eq_i]; out = set()
    for m in re.finditer(r'\(x_' + str(f) + r'\)\s*-\s*\(x_(\d+)\)', L): out.add(int(m.group(1)))
    for m in re.finditer(r'\(x_(\d+)\)\s*-\s*\(x_' + str(f) + r'\)', L): out.add(int(m.group(1)))
    return out

# start: baseline + anchors
v = base[:]
v[16742] = base[24908]; v[12186] = base[14853]
forward(v); set_quot(v)
assigned = {16742: v[16742], 12186: v[12186]}
queue = deque([16742, 12186])
cur = count(v)
print(f"start: {NEQ-len(cur)}/{NEQ} ({len(cur)} fail)")

processed = set()
rounds = 0
while queue and rounds < 4000:
    rounds += 1
    f = queue.popleft()
    fval = v[f]
    # gate partners across all eqs containing f (broken or not - align to keep checks)
    partners = set()
    for i in eqs_with(f):
        partners |= partners_in(i, f)
    for G in partners:
        if G in freeinp: continue      # only gates
        if v[G] == fval: continue       # already aligned
        w, d = additive_knob(v, G)
        if w is None: continue
        if w in assigned and assigned[w] != v[w]:
            pass
        old = v[w]
        v[w] = old + d * (fval - v[G])
        forward(v); set_quot(v)
        if v[G] != fval:
            v[w] = old; forward(v); set_quot(v); continue   # knob didn't achieve target
        new = count(v)
        if len(new) <= len(cur) + 30:   # accept (may temporarily rise; cascade fixes)
            cur = new
            if w not in assigned:
                assigned[w] = v[w]; queue.append(w)
        else:
            v[w] = old; forward(v); set_quot(v)
    if rounds % 20 == 0 or len(cur) < 12:
        print(f" round {rounds}: {NEQ-len(cur)}/{NEQ} ({len(cur)} fail) qlen={len(queue)} assigned={len(assigned)}", flush=True)
    if not cur:
        print("ALL SOLVED!"); break

cur = count(v)
print(f"FINAL: {NEQ-len(cur)}/{NEQ} ({len(cur)} fail): {sorted(cur)}")
print(f"S,T0={v[35389]%p==0},{v[6671]%p==0}  assigned={len(assigned)} free inputs")
if NEQ - len(cur) > 39016:
    json.dump({f"x_{i}": v[i] for i in range(NVARS) if v[i] != 0}, open(f'best_agentA_{NEQ-len(cur)}.json', 'w'))
    print(f"SAVED best_agentA_{NEQ-len(cur)}.json")
