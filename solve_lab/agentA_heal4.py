#!/usr/bin/env python3
"""Check control-bit config (which slack knobs are live), then aggressive heal on the
14 collateral breaks: try any free input per additive term, accept global improvement."""
import json, re, ast, sys
from agentA_harness import (p, order, gcode, definer, gates, freeinp, anc, backward_cone,
                            load_solution, forward, eqcode, eqvars, lines, NEQ, NVARS)
sys.setrecursionlimit(1000000)
boolset = set(json.load(open('boolbits.json'))['boolvars'])
base = load_solution('best/new_instance_partial_39013.json'); forward(base)

# control bit config
print("=== control config ===")
for v in [8599, 21839, 25956, 7304, 7715, 34554, 15298, 38170, 2754, 4549, 5048]:
    print(f"  x_{v} = {base[v]}  (bool={v in boolset})")
# knob liveness: x_38170=x_8599*x_21839 (gates x_21589); x_2754=x_21839*x_4549 (gates x_16787 via x_37284)
print(f"  x_16787 free? {16787 in freeinp}  gate x_2754={base[2754]} (live if !=0)")
print(f"  x_21589 free? {21589 in freeinp}  gate x_38170={base[38170]} (live if !=0)")

# Build A2 state
v = base[:]; v[16742] = base[24908]; v[12186] = base[14853]; forward(v)
def set_quot(v):
    if v[11150] % p == 0: v[30317] = -(v[11150])//p
    if (537773*v[37758]) % p == 0: v[2936] = (537773*v[37758])//p
    if v[25739] % (6672769*p) == 0: v[5146] = v[25739]//(6672769*p)
set_quot(v)

crit = {16742, 14853, 12186, 24908}
for r in [35389, 6671]:
    _, fr = backward_cone(r); crit |= fr
QUOT = {30317, 2936, 5146}
protect = set(crit) | set(QUOT)

# message-sensitive eqs for fast scoring
cd = json.load(open('agentA_code.json'))
msens = set(cd['msens'])
def count(v, subset=None):
    ns = {'__builtins__': {}, 'v': v}
    idxs = subset if subset is not None else range(NEQ)
    return set(i for i in idxs if eval(eqcode[i], ns) != 0)

# aggressive term-zeroing heal
VAR = re.compile(r'x_(\d+)')
def evn(node, v):
    if isinstance(node, ast.Constant): return node.value
    if isinstance(node, ast.Name): return v[int(node.id[2:])]
    if isinstance(node, ast.UnaryOp): return -evn(node.operand, v)
    a = evn(node.left, v); b = evn(node.right, v)
    return a+b if isinstance(node.op, ast.Add) else a-b if isinstance(node.op, ast.Sub) else a*b
def flat(node, s=1, o=None):
    if o is None: o = []
    if isinstance(node, ast.BinOp) and isinstance(node.op, (ast.Add, ast.Sub)):
        flat(node.left, s, o); flat(node.right, s*(1 if isinstance(node.op, ast.Add) else -1), o)
    else: o.append((node, s))
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
def gvars(node): return set(int(m.group(1)) for m in re.finditer(r'x_(\d+)', ast.unparse(node)))
def coeff(node, v, var):
    old = v[var]; b = evn(node, v); v[var] = old+1; c = evn(node, v)-b; v[var] = old; return c

base_fail = count(v)
print(f"\nA2 start: {NEQ-len(base_fail)}/{NEQ} ({len(base_fail)} fail)")

def heal():
    cur = count(v)
    for it in range(120):
        progressed = False
        F = sorted(cur)
        for i in F:
            root = inner(lines[i].rsplit('=', 1)[0])
            for node, s in flat(root):
                tv = evn(node, v)
                if tv == 0: continue
                cand = [x for x in gvars(node) if x in freeinp and x not in protect]
                cand.sort(key=lambda x: (x not in boolset, ))  # try non-bool slack first
                for w in cand[:25]:
                    c = coeff(node, v, w)
                    if c == 0: continue
                    if tv % c != 0: continue
                    old = v[w]; v[w] = old - tv//c
                    forward(v); set_quot(v)
                    new = count(v)
                    if len(new) < len(cur):
                        cur = new; progressed = True; break
                    v[w] = old; forward(v); set_quot(v)
                if progressed: break
            if progressed: break
        if not progressed:
            print(f" no improving move (iter {it}); stuck at {len(cur)}")
            break
        print(f" iter {it}: {NEQ-len(cur)}/{NEQ} ({len(cur)} fail)", flush=True)
        if not cur: break
    return cur

final = heal()
print(f"FINAL: {NEQ-len(final)}/{NEQ} ({len(final)} fail): {sorted(final)}")
print(f"S,T mod p zero: {v[35389]%p==0},{v[6671]%p==0}")
if NEQ-len(final) > 39016:
    json.dump({f"x_{i}": v[i] for i in range(NVARS) if v[i] != 0}, open(f'best_agentA_{NEQ-len(final)}.json', 'w'))
    print(f"SAVED best_agentA_{NEQ-len(final)}.json")
