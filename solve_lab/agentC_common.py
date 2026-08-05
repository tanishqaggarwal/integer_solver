#!/usr/bin/env python3
"""Shared setup for agentC quadratic-core solve. Copies forward_construct.py's topo forward-eval.
Exposes: p, C1, C2, gates, order, definer, gcode, forward(), val, freeinp, anc, ns, lines, eqcode,
eqvars, load_best(), CORE, and helpers to eval equation 'roots' (inner pre-square expr) mod p."""
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
val = [0]*NVARS; pinned = [False]*NVARS
for pp in A:
    vs = atom_vars(pp)
    if len(vs) == 1:
        v = next(iter(vs)); c0 = pp.get((), 0); c1 = pp.get((v,), 0); c2 = pp.get((v, v), 0)
        if c2 == 0 and c1 != 0 and (-c0) % c1 == 0 and not pinned[v]:
            val[v] = (-c0)//c1; pinned[v] = True
gate_out = set(t for t, _, _ in gates)
freeinp = set(v for v in range(NVARS) if v not in gate_out)
override = {24601: 1, 2081: 1, 30213: C2, 22162: C1, 24468: C1, 18956: C2}
for v, x in override.items():
    val[v] = x; pinned[v] = True
cand = defaultdict(list)
for gi, (t, rhs, vids) in enumerate(gates): cand[t].append(gi)
targets = set(cand); ready = [False]*NVARS
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
VAR = re.compile(r'x_(\d+)')
gcode = [compile(VAR.sub(r'v[\1]', gates[definer[order[k]]][1]), '<r>', 'eval') for k in range(len(order))]
posof = {t: k for k, t in enumerate(order)}
# free-input ancestors
anc = defaultdict(set)
for v in freeinp: anc[v] = {v}
for k, t in enumerate(order):
    s = set()
    for u in gates[definer[t]][2]: s |= anc[u]
    anc[t] = s
lines = [L for L in open('../EQUATIONS.txt').read().split('\n') if L.strip()]
eqcode = [compile(VAR.sub(r'v[\1]', L.rsplit('=', 1)[0]), '<e>', 'eval') for L in lines]
eqvars = [set(int(m) for m in VAR.findall(L)) for L in lines]
ns = {'__builtins__': {}}; ns['v'] = val
CORE = set([2071, 4573, 7123, 7469, 11854, 13660, 15299, 16622, 17726, 21382,
            22093, 25480, 25539, 28653, 29437, 31061, 32894, 32916, 34517, 34892])
# consumers for downstream partial-forward
consumers = defaultdict(list)
for k, t in enumerate(order):
    for u in gates[definer[t]][2]:
        consumers[u].append(k)

def forward():
    ns['v'] = val
    for k, t in enumerate(order): val[t] = eval(gcode[k], ns)

def downstream_ks(w):
    affected = set(); stack = [w]; seenv = set()
    while stack:
        x = stack.pop()
        if x in seenv: continue
        seenv.add(x)
        for k in consumers.get(x, ()):
            if k not in affected:
                affected.add(k); stack.append(order[k])
    return sorted(affected)

def partial_forward(ks):
    for k in ks: val[order[k]] = eval(gcode[k], ns)

def load_best(path='best/new_instance_partial_39013.json'):
    best = {int(k[2:]): v for k, v in json.load(open(path)).items()}
    for v in freeinp:
        if v in best: val[v] = best[v]
    forward()
    return best

# equation "root" = inner pre-square expression (for verifier squares E^2 the root is E)
def inner_src(lhs):
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
_rc = {}
def rootcode_of(i):
    if i not in _rc:
        _rc[i] = compile(VAR.sub(r'v[\1]', ast.unparse(inner_src(lines[i].rsplit('=', 1)[0]))), '<e>', 'eval')
    return _rc[i]

def inv(a): return pow(a % p, p - 2, p)

def is_qr(a):
    a %= p
    if a == 0: return True
    return pow(a, (p - 1)//2, p) == 1

def sqrt_mod(a):
    """Tonelli-Shanks for the field prime p; p % 4 == 3 so sqrt = a^((p+1)/4)."""
    a %= p
    if a == 0: return 0
    if not is_qr(a): return None
    r = pow(a, (p + 1)//4, p)
    assert (r*r) % p == a
    return r

if __name__ == '__main__':
    load_best()
    print(f"loaded best: NVARS={NVARS}, gates={len(gates)}, order={len(order)}, freeinp={len(freeinp)}")
    F = [i for i in range(len(lines)) if eval(eqcode[i], ns) != 0]
    print(f"satisfied {len(lines)-len(F)}/{len(lines)}; failing={len(F)}; core-fail={sorted(i for i in F if i in CORE)}")
    print(f"non-core fail={sorted(i for i in F if i not in CORE)}")
    S = val[35389] % p; T = val[6671] % p
    x3558 = val[3558] % p; x29322 = val[29322] % p; x33469 = val[33469] % p
    x27713 = val[27713] % p; x1326 = val[1326] % p; x29356 = val[29356] % p
    print(f"\n-- deep residues mod p --")
    print(f"x_3558  = {x3558}")
    print(f"x_29322 = {x29322}")
    print(f"x_29356 = {x29356}  (should be x_29322^2 = {(x29322*x29322)%p})")
    print(f"x_33469 = {x33469}")
    print(f"x_27713 = {x27713}")
    print(f"x_1326  = {x1326}")
    print(f"S = x_35389 = {S}")
    print(f"T = x_6671  = {T}")
    Scheck = (x33469*x29322*x29322 - x3558*x3558) % p
    Tcheck = (x27713*x29322 - x3558*x1326) % p
    print(f"\nS check: x_33469*x_29322^2 - x_3558^2 = {Scheck}  match={Scheck==S}")
    print(f"T check: x_27713*x_29322 - x_3558*x_1326 = {Tcheck}  match={Tcheck==T}")
    print(f"\nx_33469 is QR mod p: {is_qr(x33469)}")
    print(f"x_29322 == 0 mod p: {x29322==0}")
    print(f"regime-2 cond x_33469*x_1326^2 - x_27713^2 mod p = {(x33469*x1326*x1326 - x27713*x27713)%p}")
