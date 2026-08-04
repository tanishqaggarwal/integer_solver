#!/usr/bin/env python3
"""Accumulating exact repair from agentA_39022. Each round: repair current fails with CLEAN
(non-load-cone) handles via exact Dixon on the pivot subsystem, forward-eval, accumulate the fail
set. Track best count, save any >39022. Tests convergence vs divergence of the ripple."""
import json, re, ast, sys
from collections import defaultdict, deque
from agentE_common import build_wire, load_gates, load_lines, p, NVARS, VAR, CORE
sys.setrecursionlimit(1000000)
wire, find2, A = build_wire()
gates = load_gates()
lines = load_lines()
gate_out = set(t for t, _, _ in gates); defn = {}
for t, rhs, vids in gates:
    if t not in defn: defn[t] = tuple(vids)
freeinp = set(v for v in range(NVARS) if v not in gate_out)
anc_cache = {}
def anc(v):
    if v in anc_cache: return anc_cache[v]
    if v in freeinp or v not in defn: r = {v} if v in freeinp else set(); anc_cache[v] = r; return r
    anc_cache[v] = set(); s = set()
    for u in defn[v]: s |= anc(u)
    anc_cache[v] = s; return s
ready = [False]*NVARS
for v in range(NVARS):
    if v not in gate_out or v in freeinp: ready[v] = True
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
    if ready[t]: continue
    definer[t] = gi; order.append(t); ready[t] = True
    for gj in using[t]:
        gu[gj] -= 1
        if gu[gj] == 0: qq.append(gj)
gcode = [compile(VAR.sub(r'v[\1]', gates[definer[t]][1]), '<r>', 'eval') for t in order]
aA = {int(k[2:]): v for k, v in json.load(open('best_agentA_39022.json')).items()}
val = [0]*NVARS
for v in freeinp: val[v] = aA.get(v, 0)
ns = {'__builtins__': {}, 'v': val}
def forward():
    ns['v'] = val
    for k, t in enumerate(order): val[t] = eval(gcode[k], ns)
forward()
eqcode = [compile(VAR.sub(r'v[\1]', L.rsplit('=', 1)[0]), '<e>', 'eval') for L in lines]
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
_rc = {}
def rc(i):
    if i not in _rc: _rc[i] = compile(VAR.sub(r'v[\1]', inner_src(lines[i].rsplit('=', 1)[0])), '<r>', 'eval')
    return _rc[i]
load_cone = set()
for lv in (11150, 25739, 37758, 35389, 6671): load_cone |= anc(lv)
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
def dixon(M, b, steps=80):
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

Faccum = set()
best_cnt = 39022
F = [i for i in range(len(lines)) if eval(eqcode[i], ns)]
print(f"start: {len(lines)-len(F)} ({len(F)} fail)")
for rnd in range(12):
    Faccum |= set(F)
    FL = sorted(Faccum)
    Hset = set()
    for i in FL:
        for v in set(int(m) for m in VAR.findall(lines[i])): Hset |= anc(v)
    H = sorted(Hset - load_cone)
    base = [eval(rc(i), ns) for i in FL]
    Jac = [[0]*len(H) for _ in FL]
    for j, h in enumerate(H):
        o = val[h]; val[h] = o+1; forward()
        for ri, i in enumerate(FL): Jac[ri][j] = eval(rc(i), ns) - base[ri]
        val[h] = o
    forward()
    pr, pc = gfp_pivots(Jac); r = len(pr)
    M = [[Jac[pr[i]][pc[j]] for j in range(r)] for i in range(r)]
    rhs = [-base[pr[i]] for i in range(r)]
    y = dixon(M, rhs)
    if y is None:
        print(f"rnd {rnd}: |Facc|={len(FL)} |H|={len(H)} rank={r} -> DIXON None (inconsistent)"); break
    for j in range(r): val[H[pc[j]]] += y[j]
    forward()
    F = [i for i in range(len(lines)) if eval(eqcode[i], ns)]
    cnt = len(lines)-len(F)
    print(f"rnd {rnd}: |Facc|={len(FL)} |H|={len(H)} rank={r} -> now {cnt} ({len(F)} fail)", flush=True)
    if cnt > best_cnt:
        best_cnt = cnt
        json.dump({f"x_{i}": val[i] for i in range(NVARS)}, open(f'best_agentE_{cnt}.json', 'w'))
        print(f"  *** SAVED best_agentE_{cnt}.json ***")
    if not F: print("SOLVED!"); break
print(f"final best_cnt={best_cnt}")
