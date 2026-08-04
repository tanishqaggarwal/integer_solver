#!/usr/bin/env python3
"""Single-shot exact repair of the 26 S,T-zero breaks. Extract each break's inner linear ROOT
(strip leading const & square), build Jacobian vs union free-input cone handles, solve exactly via
mod-p pivot + Dixon p-adic lift. If Dixon returns None -> mod-p inconsistent (codeword wall).
Apply, forward-eval, count. No accumulate rounds."""
import json, re, ast, sys
from collections import defaultdict, deque
from agentE_common import build_wire, load_gates, load_lines, p, NVARS, VAR, CORE
sys.setrecursionlimit(1000000)
wire, find2, A = build_wire()
wire_set = set(wire)
gates = load_gates()
lines = load_lines()
gate_out = set(t for t, _, _ in gates); defn = {}
for t, rhs, vids in gates:
    if t not in defn: defn[t] = tuple(vids)
freeinp = set(v for v in range(NVARS) if v not in gate_out)

CORE_BRK = [11854, 29437, 32916]
NC_BRK = [3408,3841,4134,4526,5069,7276,15440,15724,15927,21600,22139,22825,27289,27999,28718,29305,31134,31269,32463,33195,36387,36390,38888]
BRK = CORE_BRK + NC_BRK

anc_cache = {}
def anc(v):
    if v in anc_cache: return anc_cache[v]
    if v in freeinp or v not in defn: r = {v} if v in freeinp else set(); anc_cache[v] = r; return r
    anc_cache[v] = set(); s = set()
    for u in defn[v]: s |= anc(u)
    anc_cache[v] = s; return s

# handles = union free-input cone (all), we will pivot to find spanning subset
Hset = set()
for i in BRK:
    for v in set(int(m) for m in VAR.findall(lines[i])): Hset |= anc(v)
Hset -= {14853, 16742}  # keep our control-zero intact
H = sorted(Hset)
print(f"{len(BRK)} breaks, {len(H)} candidate handles")

# forward-eval harness identical to stfix
PARTNERS = {30317, 5146, 2936}; CTRL = {14853, 16742}
fixed = PARTNERS | CTRL
ready = [False]*NVARS
for v in range(NVARS):
    if v not in gate_out or v in freeinp or v in fixed: ready[v] = True
gu = [0]*len(gates); using = defaultdict(list)
for gi, (t, rhs, vids) in enumerate(gates):
    u = 0
    for v in vids:
        if not ready[v]: u += 1
        using[v].append(gi)
    gu[gi] = u
definer = {}; order = []
qq = deque(gi for gi in range(len(gates)) if gu[gi] == 0)
while qq:
    gi = qq.popleft(); t, rhs, vids = gates[gi]
    if ready[t] or t in fixed: continue
    definer[t] = gi; order.append(t); ready[t] = True
    for gj in using[t]:
        gu[gj] -= 1
        if gu[gj] == 0: qq.append(gj)
gcode = [compile(VAR.sub(r'v[\1]', gates[definer[t]][1]), '<r>', 'eval') for t in order]
best = {int(k[2:]): v for k, v in json.load(open('best/new_instance_partial_39013.json')).items()}
val = [0]*NVARS
for v in freeinp: val[v] = best.get(v, 0)
ns = {'__builtins__': {}, 'v': val}
def forward():
    ns['v'] = val
    for k, t in enumerate(order): val[t] = eval(gcode[k], ns)
forward()
val[14853] = val[12186]; val[16742] = val[24908]; forward()
eqcode = [compile(VAR.sub(r'v[\1]', L.rsplit('=', 1)[0]), '<e>', 'eval') for L in lines]

# inner-root extraction (strip const* and square) -> compiled code for the linear form
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
    return ast.unparse(node)
rootcode = [compile(VAR.sub(r'v[\1]', inner_src(lines[i].rsplit('=', 1)[0])), '<r>', 'eval') for i in BRK]
def roots(): ns['v'] = val; return [eval(c, ns) for c in rootcode]

def inv(a): return pow(a % p, p-2, p)
def matinv(M):
    r = len(M); Aug = [[M[i][j] % p for j in range(r)] + [1 if j == i else 0 for j in range(r)] for i in range(r)]
    for c in range(r):
        piv = next((i for i in range(c, r) if Aug[i][c] % p != 0), None)
        if piv is None: return None
        Aug[c], Aug[piv] = Aug[piv], Aug[c]; iv = inv(Aug[c][c]); Aug[c] = [(x*iv) % p for x in Aug[c]]
        for i in range(r):
            if i != c and Aug[i][c] % p != 0:
                f = Aug[i][c]; Aug[i] = [(Aug[i][k]-f*Aug[c][k]) % p for k in range(2*r)]
    return [[Aug[i][r+j] for j in range(r)] for i in range(r)]
def dixon(M, b, steps=40):
    r = len(M); Mi = matinv(M)
    if Mi is None: return None
    x = [0]*r; bb = b[:]; mod = 1
    for _ in range(steps):
        bm = [bb[i] % p for i in range(r)]
        xi = [sum(Mi[i][k]*bm[k] for k in range(r)) % p for i in range(r)]
        for i in range(r): x[i] += mod*xi[i]
        nb = []
        for i in range(r):
            s = bb[i] - sum(M[i][k]*xi[k] for k in range(r))
            if s % p != 0: return None
            nb.append(s//p)
        bb = nb; mod *= p
        if all(z == 0 for z in bb): break
    half = mod//2; y = []
    for xi in x:
        xi %= mod
        if xi > half: xi -= mod
        y.append(xi)
    return y

base = roots()
print(f"nonzero roots at start: {sum(1 for r in base if r)}/{len(base)}")
# Jacobian of roots vs handles
Jac = [[0]*len(H) for _ in BRK]
for j, h in enumerate(H):
    o = val[h]; val[h] = o+1; forward(); r1 = roots()
    for ri in range(len(BRK)): Jac[ri][j] = r1[ri] - base[ri]
    val[h] = o
forward()
# mod-p pivots to find spanning columns/rows
def gfp_pivots(M):
    m = len(M); n = len(M[0]) if m else 0
    Mx = [[M[i][j] % p for j in range(n)] for i in range(m)]
    rowmap = list(range(m)); pr = 0; pivr = []; pivc = []
    for c in range(n):
        piv = next((i for i in range(pr, m) if Mx[i][c] % p != 0), None)
        if piv is None: continue
        Mx[pr], Mx[piv] = Mx[piv], Mx[pr]; rowmap[pr], rowmap[piv] = rowmap[piv], rowmap[pr]
        iv = inv(Mx[pr][c]); Mx[pr] = [(x*iv) % p for x in Mx[pr]]
        for i in range(m):
            if i != pr and Mx[i][c] % p != 0:
                f = Mx[i][c]; Mx[i] = [(Mx[i][k]-f*Mx[pr][k]) % p for k in range(n)]
        pivr.append(rowmap[pr]); pivc.append(c); pr += 1
        if pr >= m: break
    return pivr, pivc
pr, pc = gfp_pivots(Jac)
r = len(pr)
print(f"Jacobian rank (mod p) = {r} of {len(BRK)} breaks")
M = [[Jac[pr[i]][pc[j]] for j in range(r)] for i in range(r)]
rhs = [-base[pr[i]] for i in range(r)]
y = dixon(M, rhs)
if y is None:
    print("DIXON FAILED: mod-p inconsistent (codeword wall) on the pivot subsystem")
else:
    for j in range(r): val[H[pc[j]]] += y[j]
    forward()
    F = [i for i in range(len(lines)) if eval(eqcode[i], ns)]
    c = [i for i in F if i in CORE]; nc = [i for i in F if i not in CORE]
    print(f"AFTER single-shot repair: {len(lines)-len(F)}/{len(lines)} ({len(F)} fail); core {len(c)}, noncore {len(nc)}")
    print(f"  remaining fails: {sorted(F)[:30]}")
    if len(lines)-len(F) > 39013:
        json.dump({f"x_{i}": val[i] for i in range(NVARS)}, open(f'best_agentE_{len(lines)-len(F)}.json', 'w'))
        print(f"  SAVED best_agentE_{len(lines)-len(F)}.json")
