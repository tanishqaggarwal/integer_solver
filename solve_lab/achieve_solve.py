#!/usr/bin/env python3
"""Solve by recursive circuit inversion.

To make an output var equal a target, recurse through its gate to the free inputs:
  copy x=a          -> achieve(a, target)
  sum  x=a+b        -> achieve(a, target), achieve(b, 0)
  diff x=a-b        -> achieve(a, target), achieve(b, 0)
  prod x=a*b, t!=0  -> put target on the free-input side, other factor = 1
  prod x=a*b, t==0  -> zero the easier factor
  const c           -> require c==target
Free inputs get set. Then forward-evaluate the whole circuit and verify.
Targets: x_37892 = BIGCONST, x_13682 = H2 (routes the gaps through product slacks)."""
import json, re, ast, sys
from collections import defaultdict, deque
from propagate import load_atoms, atom_vars, NVARS

A = load_atoms()
base = {int(k[2:]): x for k, x in json.load(open('rebuilt_partial.json')).items()}
V0 = base[23917]; H1 = abs(A[602].get((), 0)); H2 = abs(A[1465].get((), 0)); BIGCONST = H1 // 8863713
gates = {}
gates_all = []
with open('atoms/gates.jsonl') as f:
    for line in f:
        d = json.loads(line); gates_all.append((d['t'], d['rhs'], tuple(d['vids'])))
        if d['t'] not in gates: gates[d['t']] = ast.parse(d['rhs'], mode='eval').body
summ = json.load(open('atoms/summary.json')); inputs = set(summ['inputs'])
# also treat pinned single-var-atom vars as fixed
pin = {}
for p in A:
    vs = atom_vars(p)
    if len(vs) == 1:
        v = next(iter(vs)); c0 = p.get((), 0); c1 = p.get((v,), 0); c2 = p.get((v, v), 0)
        if c2 == 0 and c1 != 0 and (-c0) % c1 == 0 and v not in pin:
            pin[v] = (-c0) // c1

setv = {}      # free-input assignments chosen by inversion
conflicts = []

# does this subtree contain a free input reachable via copy/sum/diff (an "easy" side)?
def free_reachable(node, depth=0):
    if depth > 40: return False
    if isinstance(node, ast.Name):
        v = int(node.id[2:])
        if v in inputs: return True
        if v in gates: return free_reachable(gates[v], depth+1)
        return False
    if isinstance(node, ast.BinOp):
        return free_reachable(node.left, depth+1) or free_reachable(node.right, depth+1)
    return False

def achieve(node, target, depth=0):
    if depth > 60:
        conflicts.append(('depth', target)); return
    if isinstance(node, ast.Constant):
        if node.value != target: conflicts.append(('const', node.value, target))
        return
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
        achieve(node.operand, -target, depth+1); return
    if isinstance(node, ast.Name):
        v = int(node.id[2:])
        if v in pin:
            if pin[v] != target: conflicts.append(('pin', v, pin[v], target))
            return
        if v in inputs or v not in gates:
            if v in setv and setv[v] != target: conflicts.append(('inp', v, setv[v], target))
            else: setv[v] = target
            return
        achieve(gates[v], target, depth+1); return
    if isinstance(node, ast.BinOp):
        a, b = node.left, node.right
        if isinstance(node.op, ast.Add):
            achieve(a, target, depth+1); achieve(b, 0, depth+1); return
        if isinstance(node.op, ast.Sub):
            achieve(a, target, depth+1); achieve(b, 0, depth+1); return
        if isinstance(node.op, ast.Mult):
            if target == 0:
                # zero the side that's easier to zero (free reachable), default left
                if free_reachable(a): achieve(a, 0, depth+1)
                else: achieve(b, 0, depth+1)
            else:
                # put target on the free-reachable side, 1 on the other
                if free_reachable(b) and not free_reachable(a):
                    achieve(b, target, depth+1); achieve(a, 1, depth+1)
                else:
                    achieve(a, target, depth+1); achieve(b, 1, depth+1)
            return
    conflicts.append(('unhandled', ast.dump(node)[:40]))

# route the two gaps
achieve(gates[37892], BIGCONST)
achieve(gates[13682], H2)
print(f"inversion set {len(setv)} free inputs; conflicts: {len(conflicts)}: {conflicts[:6]}", flush=True)

# forward-eval everything from base + the chosen inputs
val = [0]*NVARS
for v, x in pin.items(): val[v] = x
for v, x in setv.items(): val[v] = x
readyset = set(pin) | set(setv)
cand = defaultdict(list)
for gi,(t,rhs,vids) in enumerate(gates_all): cand[t].append(gi)
targets = set(cand)
ready = [False]*NVARS
for v in range(NVARS):
    if v not in targets or v in readyset: ready[v] = True
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
for t in order:
    ns['v'] = val; val[t] = eval(code[t], ns)

json.dump({f"x_{i}": val[i] for i in range(NVARS)}, open('achieved.json','w'))
def ev(poly):
    s=0
    for m,c in poly.items():
        t=c
        for x in m: t*=val[x]
        s+=t
    return s
nz=[i for i in range(len(A)) if ev(A[i])!=0]
print(f"nonzero atoms: {len(nz)} -> {nz[:24]}", flush=True)
