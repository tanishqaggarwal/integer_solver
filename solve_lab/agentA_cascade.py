#!/usr/bin/env python3
"""Cascade alignment healer. For each nonzero additive term in a failing eq, candidate
knobs = direct free vars + live additive free knobs of gates in the term (recursively
found). Set knob to zero the term; accept if global fail count does not increase."""
import json, re, ast, sys
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

# --- additive-knob cache: for a gate, a free input w with d(gate)/dw = +-1 ---
_knob = {}
def free_cone(root, limit=400):
    seen=set(); st=[root]
    while st and len(seen)<limit:
        u=st.pop()
        if u in seen: continue
        seen.add(u)
        for w in gvids.get(u, ()):
            if w not in seen: st.append(w)
    return [u for u in seen if u in freeinp]
def additive_knob(v, gate, avoid):
    if gate in _knob:
        w, d = _knob[gate]
        if w is None or w in avoid: return None, None
        return w, d
    g0 = v[gate]
    res = (None, None)
    for w in free_cone(gate):
        if w in protect: continue
        old = v[w]; v[w] = old+1; forward(v); d = v[gate]-g0; v[w] = old
        if d in (1, -1):
            res = (w, d); break
    forward(v)
    _knob[gate] = res
    w, d = res
    if w is None or w in avoid: return None, None
    return w, d

# --- term flattening ---
def evn(node, v):
    if isinstance(node, ast.Constant): return node.value
    if isinstance(node, ast.Name): return v[int(node.id[2:])]
    if isinstance(node, ast.UnaryOp): return -evn(node.operand, v)
    a = evn(node.left, v); b = evn(node.right, v)
    return a+b if isinstance(node.op, ast.Add) else a-b if isinstance(node.op, ast.Sub) else a*b
def flatterms(node, o=None):
    if o is None: o = []
    if isinstance(node, ast.BinOp) and isinstance(node.op, (ast.Add, ast.Sub)):
        flatterms(node.left, o); flatterms(node.right, o)
    else: o.append(node)
    return o
def inner(lhs):
    node = ast.parse(lhs, mode='eval').body
    while isinstance(node, ast.BinOp) and isinstance(node.op, ast.Mult):
        a, b = node.left, node.right
        ca = isinstance(a, ast.Constant) or (isinstance(a, ast.UnaryOp) and isinstance(a.operand, ast.Constant))
        cb = isinstance(b, ast.Constant) or (isinstance(b, ast.UnaryOp) and isinstance(b.operand, ast.Constant))
        if ca and not cb: node = b
        elif cb and not ca: node = a
        elif ast.unparse(a) == ast.unparse(b): node = a
        else: break
    return node
_rootcache = {}
def rootterms(i):
    if i not in _rootcache: _rootcache[i] = flatterms(inner(lines[i].rsplit('=', 1)[0]))
    return _rootcache[i]
def gvars(node): return set(int(m.group(1)) for m in re.finditer(r'x_(\d+)', ast.unparse(node)))

# --- build start state: baseline + core anchors ---
v = base[:]
v[16742] = base[24908]; v[12186] = base[14853]
forward(v); set_quot(v)
cur = count(v)
print(f"start (anchors): {NEQ-len(cur)}/{NEQ} ({len(cur)} fail)")

def try_zero_term(node, v):
    """try to make additive term `node` == 0 via a knob (direct free or gate additive knob)."""
    tv = evn(node, v)
    if tv == 0: return False
    # candidate knobs: direct free vars in term
    cand = []
    for x in gvars(node):
        if x in freeinp and x not in protect: cand.append((x, None))  # None => compute coeff
    # gate additive knobs
    for x in gvars(node):
        if x not in freeinp:
            w, d = additive_knob(v, x, protect)
            if w is not None: cand.append((w, ('gate', x, d)))
    for w, meta in cand:
        old = v[w]
        # compute coefficient of w in the term value
        v[w] = old + 1; tv1 = evn(node, v); v[w] = old
        c = tv1 - tv
        if c == 0: continue
        if tv % c != 0: continue
        v[w] = old - tv // c
        forward(v); set_quot(v)
        new = count(v)
        if len(new) <= len(cur_holder[0]):
            cur_holder[0] = new
            return True
        v[w] = old; forward(v); set_quot(v)
    return False

cur_holder = [cur]
for it in range(300):
    progressed = False
    for i in sorted(cur_holder[0]):
        for node in rootterms(i):
            before = len(cur_holder[0])
            if try_zero_term(node, v):
                if len(cur_holder[0]) < before:
                    progressed = True; break
        if progressed: break
    c = cur_holder[0]
    if it % 1 == 0:
        print(f" iter {it}: {NEQ-len(c)}/{NEQ} ({len(c)} fail)  S,T0={v[35389]%p==0},{v[6671]%p==0}", flush=True)
    if not c:
        print("ALL SOLVED!"); break
    if not progressed:
        print(f" stuck at {len(c)}: {sorted(c)}"); break
c = count(v)
print(f"FINAL: {NEQ-len(c)}/{NEQ} ({len(c)} fail): {sorted(c)}")
if NEQ - len(c) > 39016:
    json.dump({f"x_{i}": v[i] for i in range(NVARS) if v[i] != 0}, open(f'best_agentA_{NEQ-len(c)}.json', 'w'))
    print(f"SAVED best_agentA_{NEQ-len(c)}.json")
