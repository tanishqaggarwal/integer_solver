#!/usr/bin/env python3
"""Healer: seed with Route A (S,T=0), protect S,T-critical free inputs, run
forward_construct-style try_set on the failing equations to heal wiring damage."""
import json, re, ast, sys
from agentA_harness import (p, order, gcode, definer, gates, freeinp, anc, backward_cone,
                            load_solution, forward, eqcode, eqvars, lines, NEQ, NVARS)
sys.setrecursionlimit(1000000)

base = load_solution('best/new_instance_partial_39013.json'); forward(base)
val = base[:]
# Route A pins
val[16742] = base[24908]; val[14853] = base[12186]
forward(val)

# protect S,T-critical free inputs (so healing can't un-zero S,T)
crit = {16742, 14853, 12186}
for r in [35389, 6671]:
    _, fr = backward_cone(r); crit |= fr
QUOT = {30317, 2936, 5146}

ns = {'__builtins__': {}}
def set_quot():
    if val[11150] % p == 0: val[30317] = -(val[11150])//p
    if (537773*val[37758]) % p == 0: val[2936] = (537773*val[37758])//p
    if val[25739] % (6672769*p) == 0: val[5146] = val[25739]//(6672769*p)

# ---- ported try_set machinery from forward_construct ----
def evn(node):
    if isinstance(node, ast.Constant): return node.value
    if isinstance(node, ast.Name): return val[int(node.id[2:])]
    if isinstance(node, ast.UnaryOp): return -evn(node.operand)
    a = evn(node.left); b = evn(node.right)
    return a+b if isinstance(node.op, ast.Add) else a-b if isinstance(node.op, ast.Sub) else a*b
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
def flat(node, s=1, o=None):
    if o is None: o = []
    if isinstance(node, ast.BinOp) and isinstance(node.op, (ast.Add, ast.Sub)):
        flat(node.left, s, o); flat(node.right, s*(1 if isinstance(node.op, ast.Add) else -1), o)
    else: o.append(node)
    return o
astcache = {}
def rootast(i):
    if i not in astcache: astcache[i] = inner(lines[i].rsplit('=', 1)[0])
    return astcache[i]
def gvars(node): return set(int(m.group(1)) for m in re.finditer(r'x_(\d+)', ast.unparse(node)))
def free_deps(node):
    s = set()
    for v in gvars(node): s |= anc.get(v, {v} if v in freeinp else set())
    return s
def coeff(node, v):
    old = val[v]; b = evn(node); val[v] = old+1; c = evn(node)-b; val[v] = old; return c

determined = set(crit) | set(QUOT)
def dep_final(w): return w in determined or val[w] == 0
def try_set(t):
    frees = [v for v in gvars(t) if v in freeinp and v not in determined]
    if not frees: return False
    rem = None; quot = None
    others_det = lambda v: all(dep_final(w) for w in free_deps(t)-{v})
    for v in frees:
        if not others_det(v): continue
        c = coeff(t, v)
        if c == 0: continue
        if c % p == 0:
            if quot is None: quot = (v, c)
        else:
            if rem is None: rem = (v, c)
    cur = evn(t)
    if cur == 0: return False
    if rem is not None and quot is not None:
        vr, cr = rem; vq, cq = quot
        if all(dep_final(w) for w in free_deps(t)-{vr, vq}):
            dr = (-cur*pow(cr % p, p-2, p)) % p
            Rp = cur + cr*dr
            if Rp % cq == 0:
                val[vr] += dr; val[vq] += (-Rp)//cq; determined.add(vr); determined.add(vq); return True
    if rem is not None:
        vr, cr = rem
        if cur % cr == 0: val[vr] -= cur//cr; determined.add(vr); return True
    if quot is not None:
        vq, cq = quot
        if cur % cq == 0: val[vq] -= cur//cq; determined.add(vq); return True
    return False

set_quot(); forward(val); set_quot()
ns['v'] = val
F0 = [i for i in range(NEQ) if eval(eqcode[i], ns) != 0]
print(f"start heal: {NEQ-len(F0)}/{NEQ} ({len(F0)} fail)")
for it in range(200):
    ns['v'] = val
    F = [i for i in range(NEQ) if eval(eqcode[i], ns) != 0]
    if it % 10 == 0 or len(F) < 25:
        print(f" iter {it}: {NEQ-len(F)}/{NEQ} ({len(F)} fail) det={len(determined)}", flush=True)
    if not F: print("ALL SOLVED!"); break
    changed = False
    for i in F:
        for t in flat(rootast(i)):
            if evn(t) != 0 and try_set(t): changed = True
    forward(val); set_quot()
    if not changed:
        print(f" stuck at {len(F)} fail: {sorted(F)[:30]}"); break
ns['v'] = val
F = [i for i in range(NEQ) if eval(eqcode[i], ns) != 0]
print(f"FINAL heal: {NEQ-len(F)}/{NEQ} ({len(F)} fail): {sorted(F)[:30]}")
# verify S,T still 0
print(f"S mod p={val[35389]%p==0}  T mod p={val[6671]%p==0}")
if len(F) < 20:
    out = {f"x_{i}": val[i] for i in range(NVARS) if val[i] != 0}
    json.dump(out, open(f'best_agentA_heal_{NEQ-len(F)}.json', 'w'))
    print(f"SAVED best_agentA_heal_{NEQ-len(F)}.json")
