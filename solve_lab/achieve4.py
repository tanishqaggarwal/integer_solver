#!/usr/bin/env python3
"""Recursive circuit inversion with target-tracking coordination.

req(v,t): require var v == t. Records target[v] (conflict if inconsistent), then solves
v's gate expression = t. Key rules:
 - Mult target!=0: route target to a factor that can carry it (a free input reachable via
   value ops), set the other factor to 1.
 - Mult target==0: zero the factor NOT already targeted nonzero.
 - Add/Sub: route target to the value-carrying side, 0 to the other.
Coordination via target[] makes OR-gates (a+b-a*b) resolve consistently."""
import json, re, ast, sys
from collections import defaultdict, deque
from propagate import load_atoms, atom_vars, NVARS

A = load_atoms()
base = {int(k[2:]): x for k, x in json.load(open('rebuilt_partial.json')).items()}
H1 = abs(A[602].get((), 0)); H2 = abs(A[1465].get((), 0)); BIGCONST = H1 // 8863713
gates = {}; gates_all = []
with open('atoms/gates.jsonl') as f:
    for line in f:
        d = json.loads(line); gates_all.append((d['t'], d['rhs'], tuple(d['vids'])))
        if d['t'] not in gates: gates[d['t']] = ast.parse(d['rhs'], mode='eval').body
summ = json.load(open('atoms/summary.json')); inputs = set(summ['inputs'])
pin = {}
for p in A:
    vs = atom_vars(p)
    if len(vs) == 1:
        v = next(iter(vs)); c0 = p.get((), 0); c1 = p.get((v,), 0); c2 = p.get((v, v), 0)
        if c2 == 0 and c1 != 0 and (-c0) % c1 == 0 and v not in pin:
            pin[v] = (-c0) // c1

target = {}; conflicts = []
sys.setrecursionlimit(100000)

def name_var(node):
    return int(node.id[2:]) if isinstance(node, ast.Name) else None

# can this node carry a nonzero value down to a free input (value path: copy/sum/sub, or
# one factor of a product)? returns True if a free input is reachable to hold magnitude.
def carries(node, depth=0):
    if depth > 50: return False
    if isinstance(node, ast.Constant): return False
    if isinstance(node, ast.UnaryOp): return carries(node.operand, depth+1)
    if isinstance(node, ast.Name):
        v = int(node.id[2:])
        if v in pin: return False
        if v in inputs or v not in gates: return True
        return carries(gates[v], depth+1)
    if isinstance(node, ast.BinOp):
        if isinstance(node.op, ast.Mult):
            return carries(node.left, depth+1) or carries(node.right, depth+1)
        return carries(node.left, depth+1) or carries(node.right, depth+1)
    return False

def req(v, t):
    if v in pin:
        if pin[v] != t: conflicts.append(('pin', v, pin[v], t))
        return
    if v in target:
        if target[v] != t: conflicts.append(('var', v, target[v], t))
        return
    target[v] = t
    if v in inputs or v not in gates:
        return
    solve(gates[v], t)

def solve(node, t, depth=0):
    if depth > 400: conflicts.append(('depth',)); return
    if isinstance(node, ast.Constant):
        if node.value != t: conflicts.append(('const', node.value, t)); return
        return
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
        solve(node.operand, -t, depth+1); return
    if isinstance(node, ast.Name):
        req(int(node.id[2:]), t); return
    if isinstance(node, ast.BinOp):
        a, b = node.left, node.right
        if isinstance(node.op, (ast.Add, ast.Sub)):
            sgn = 1 if isinstance(node.op, ast.Add) else -1
            # route t to the side that carries; 0 to the other
            if carries(a): solve(a, t, depth+1); solve(b, 0, depth+1)
            elif carries(b): solve(b, sgn*t, depth+1); solve(a, 0, depth+1)
            else: solve(a, t, depth+1); solve(b, 0, depth+1)
            return
        if isinstance(node.op, ast.Mult):
            va, vb = name_var(a), name_var(b)
            if t == 0:
                # zero the factor not already targeted nonzero
                if va is not None and target.get(va, 0) != 0: solve(b, 0, depth+1)
                elif vb is not None and target.get(vb, 0) != 0: solve(a, 0, depth+1)
                elif carries(a): solve(a, 0, depth+1)
                else: solve(b, 0, depth+1)
            else:
                # put t on the carrying side, 1 on the other
                if carries(a) and not carries(b): solve(a, t, depth+1); solve(b, 1, depth+1)
                elif carries(b): solve(b, t, depth+1); solve(a, 1, depth+1)
                else: solve(a, t, depth+1); solve(b, 1, depth+1)
            return
    conflicts.append(('unhandled',))

# Phase 1: activate x_7715=1 via a CLEAN load bit (its load absorber avoids the verifier squares)
target[36314]=1
print(f"activation: targets {len(target)}; free inputs {sum(1 for v in target if v in inputs or v not in gates)}; conflicts {len(conflicts)}: {conflicts[:6]}", flush=True)
# Phase 2: with x_15298=1 the slack products pass values through; set the value inputs and gap outputs
for v, x in [(16742, BIGCONST), (12186, H2), (18956, BIGCONST), (24468, H2)]:
    if v not in target: target[v] = x

# forward-eval from base(=0 free) overridden by chosen free-input targets
val = [0]*NVARS
for v, x in pin.items(): val[v] = x
setfree = {v: target[v] for v in target if (v in inputs or v not in gates) and v not in pin}
for v, x in setfree.items(): val[v] = x
readyset = set(pin) | set(setfree)
cand = defaultdict(list)
for gi,(t,rhs,vids) in enumerate(gates_all): cand[t].append(gi)
targets_set = set(cand)
ready = [False]*NVARS
for v in range(NVARS):
    if v not in targets_set or v in readyset: ready[v] = True
gu = [0]*len(gates_all); using = defaultdict(list)
for gi,(t,rhs,vids) in enumerate(gates_all):
    u = 0
    for v in vids:
        if not ready[v]: u += 1
        using[v].append(gi)
    gu[gi] = u
definer = {}; order = []
q = deque(gi for gi in range(len(gates_all)) if gu[gi]==0)
while q:
    gi = q.popleft(); t,rhs,vids = gates_all[gi]
    if ready[t]: continue
    definer[t] = gi; order.append(t); ready[t] = True
    for gj in using[t]:
        gu[gj] -= 1
        if gu[gj]==0: q.append(gj)
VAR = re.compile(r'x_(\d+)')
code = {t: compile(VAR.sub(r'v[\1]', gates_all[definer[t]][1]), '<r>', 'eval') for t in order}
ns = {'__builtins__': {}}
def ev(poly):
    s=0
    for m,c in poly.items():
        z=c
        for x in m: z*=val[x]
        s+=z
    return s
avar = [atom_vars(A[i]) for i in range(len(A))]
# iterative forward-eval + load-absorption: each broken load atom is absorbed by setting a
# free-input variable (the load absorber) to satisfy it.
freeset = set(v for v in range(NVARS) if v not in targets_set) - set(setfree) - set(pin)  # absorbers only
for it in range(60):
    for t in order:
        ns['v'] = val; val[t] = eval(code[t], ns)
    changed = False
    for ai in range(len(A)):
        r = ev(A[ai])
        if r == 0: continue
        for v in avar[ai]:
            if v not in freeset or v in pin: continue
            d = 0; lin = True
            for m, c in A[ai].items():
                if v in m:
                    if m.count(v) >= 2: lin = False; break
                    z = c
                    for x in m:
                        if x != v: z *= val[x]
                    d += z
            if lin and d != 0 and r % d == 0:
                val[v] -= r // d; changed = True; break
    if not changed: break
nz=[i for i in range(len(A)) if ev(A[i])!=0]
json.dump({f"x_{i}": val[i] for i in range(NVARS)}, open('achieved4.json','w'))
print(f"after {it+1} repair rounds: nonzero atoms: {len(nz)} -> {nz[:24]}", flush=True)
