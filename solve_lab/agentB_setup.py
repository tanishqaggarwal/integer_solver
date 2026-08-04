#!/usr/bin/env python3
"""Shared setup for AgentB global mod-p solve.
Builds: topological gate order, free-input set, poly form of every gate rhs and every
equation-root (via inner()), mod-p forward-eval, and exact tangent-linear (symbolic
gradient) propagation mod p.  Heavy parsing is cached to scratchpad.
"""
import json, re, ast, sys, os, pickle, time
from collections import defaultdict, deque

sys.setrecursionlimit(1000000)
p = 2**256 - 2**32 - 977
HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = '/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/agentB_cache.pkl'
NVARS = 38748

hc = json.load(open(HERE + '/huge_consts.json')); C1 = int(hc['C1']); C2 = int(hc['C2'])
override = {24601: 1, 2081: 1, 30213: C2, 22162: C1, 24468: C1, 18956: C2}
DEFAULT_OVERRIDE = dict(override)

VAR = re.compile(r'x_(\d+)')

def ast_to_poly(node):
    """Return dict{ sorted-tuple-of-varids : int coef } for an AST expression node."""
    if isinstance(node, ast.Expression):
        return ast_to_poly(node.body)
    if isinstance(node, ast.Constant):
        return {(): int(node.value)} if node.value != 0 else {}
    if isinstance(node, ast.Name):
        return {(int(node.id[2:]),): 1}
    if isinstance(node, ast.UnaryOp):
        if isinstance(node.op, ast.USub):
            return {m: -c for m, c in ast_to_poly(node.operand).items()}
        if isinstance(node.op, ast.UAdd):
            return ast_to_poly(node.operand)
        raise ValueError('unary')
    if isinstance(node, ast.BinOp):
        a = ast_to_poly(node.left); b = ast_to_poly(node.right)
        if isinstance(node.op, ast.Add):
            out = defaultdict(int, a)
            for m, c in b.items(): out[m] += c
            return {m: c for m, c in out.items() if c != 0}
        if isinstance(node.op, ast.Sub):
            out = defaultdict(int, a)
            for m, c in b.items(): out[m] -= c
            return {m: c for m, c in out.items() if c != 0}
        if isinstance(node.op, ast.Mult):
            out = defaultdict(int)
            for m1, c1 in a.items():
                for m2, c2 in b.items():
                    out[tuple(sorted(m1 + m2))] += c1 * c2
            return {m: c for m, c in out.items() if c != 0}
        raise ValueError('binop ' + str(node.op))
    raise ValueError('node ' + str(node))

def inner(lhs):
    """Strip outer constant-multipliers and squares to get the equation ROOT node."""
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

def build(ov=None, cache=None):
    if ov is None: ov = DEFAULT_OVERRIDE
    if cache is None: cache = CACHE
    t0 = time.time()
    # atoms single-var pins
    from propagate import load_atoms, atom_vars
    A = load_atoms()
    pins = {}
    for pp in A:
        vs = atom_vars(pp)
        if len(vs) == 1:
            v = next(iter(vs)); c0 = pp.get((), 0); c1 = pp.get((v,), 0); c2 = pp.get((v, v), 0)
            if c2 == 0 and c1 != 0 and (-c0) % c1 == 0 and v not in pins:
                pins[v] = (-c0) // c1
    gates = []
    with open(HERE + '/atoms/gates.jsonl') as f:
        for line in f:
            d = json.loads(line); gates.append((d['t'], d['rhs'], tuple(d['vids'])))
    gate_out = set(t for t, _, _ in gates)
    freeinp = set(v for v in range(NVARS) if v not in gate_out)
    # topo order (replicate forward_construct)
    pinned = set(pins) | set(ov)
    cand = defaultdict(list)
    for gi, (t, rhs, vids) in enumerate(gates): cand[t].append(gi)
    targets = set(cand); ready = [False] * NVARS
    for v in range(NVARS):
        if v not in targets or v in freeinp or v in pinned: ready[v] = True
    gu = [0] * len(gates); using = defaultdict(list)
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
    # parse gate polys in topo order
    gate_poly = []  # list of (t, polydict) in topo order
    for t in order:
        _, rhs, _ = gates[definer[t]]
        pol = ast_to_poly(ast.parse(rhs, mode='eval'))
        gate_poly.append((t, pol))
    # equations -> root polys
    lines = [L for L in open(HERE + '/../EQUATIONS.txt').read().split('\n') if L.strip()]
    root_poly = []
    eqvars = []
    for L in lines:
        lhs = L.rsplit('=', 1)[0]
        rp = ast_to_poly(inner(lhs))
        root_poly.append(rp)
        eqvars.append(frozenset(int(m) for m in VAR.findall(lhs)))
    data = dict(pins=pins, freeinp=sorted(freeinp), order=order,
                gate_poly=gate_poly, root_poly=root_poly, eqvars=eqvars,
                gate_out=sorted(gate_out), override=dict(ov))
    with open(cache, 'wb') as f:
        pickle.dump(data, f)
    print(f"[setup] built+cached ({cache}) in {time.time()-t0:.1f}s", file=sys.stderr)
    return data

def load(ov=None, cache=None):
    if ov is None: ov = DEFAULT_OVERRIDE
    if cache is None: cache = CACHE
    if os.path.exists(cache):
        with open(cache, 'rb') as f:
            return pickle.load(f)
    return build(ov, cache)

# ---- runtime helpers ----
class Env:
    def __init__(self, data):
        self.freeinp = data['freeinp']
        self.freeset = set(self.freeinp)
        self.order = data['order']
        self.gate_poly = data['gate_poly']
        self.root_poly = data['root_poly']
        self.eqvars = data['eqvars']
        self.pins = data['pins']
        self.override = data.get('override', DEFAULT_OVERRIDE)
        self.gate_out = set(data['gate_out'])
        self.valp = [0] * NVARS  # mod p values
        self.forced = {}         # var -> value (treated as a mod-p CONSTANT, gradient 0)

    def set_from_solution(self, sol_int):
        """sol_int: dict var->int (full). Set free inputs' valp; pins/override applied."""
        vp = self.valp
        for v in range(NVARS):
            vp[v] = 0
        for v, x in self.pins.items(): vp[v] = x % p
        for v, x in self.override.items(): vp[v] = x % p
        for v in self.freeset:
            if v in sol_int: vp[v] = sol_int[v] % p
        # pins/override may be free inputs; ensure applied last for free inputs too
        for v, x in self.override.items(): vp[v] = x % p
        self.forward()

    def set_free(self, freevals):
        """freevals: dict var->int for free inputs only (others left/pins)."""
        vp = self.valp
        for v in range(NVARS): vp[v] = 0
        for v, x in self.pins.items(): vp[v] = x % p
        for v in self.freeset:
            vp[v] = freevals.get(v, 0) % p
        for v, x in self.override.items(): vp[v] = x % p
        self.forward()

    def forward(self):
        vp = self.valp
        forced = self.forced
        if forced:
            for v, x in forced.items(): vp[v] = x % p
        for t, pol in self.gate_poly:
            if forced and t in forced: continue
            s = 0
            for m, c in pol.items():
                term = c
                for v in m: term *= vp[v]
                s += term
            vp[t] = s % p

    def root_val(self, i):
        vp = self.valp; s = 0
        for m, c in self.root_poly[i].items():
            term = c
            for v in m: term *= vp[v]
            s += term
        return s % p

    def all_root_residuals(self):
        """Return dict i->residual (mod p) for nonzero roots."""
        out = {}
        vp = self.valp
        for i, rp in enumerate(self.root_poly):
            s = 0
            for m, c in rp.items():
                term = c
                for v in m: term *= vp[v]
                s += term
            s %= p
            if s: out[i] = s
        return out

    def tangent_linear(self):
        """Compute grad[v] = {freeinput: coef} mod p for every var, at current valp.
        Exact mod-p Jacobian of forward-eval wrt free inputs."""
        vp = self.valp
        grad = [None] * NVARS
        jac_free = getattr(self, 'jac_free', set())
        for v in self.freeset:
            grad[v] = {v: 1}
        for v in self.forced:
            grad[v] = {v: 1} if v in jac_free else {}   # forced-but-solvable -> column; else constant
        for t, pol in self.gate_poly:
            if t in self.forced:
                grad[t] = {t: 1} if t in jac_free else {}
                continue
            g = {}
            self._accum_poly_grad(pol, g, grad, vp)
            grad[t] = {k: c for k, c in g.items() if c}
        self.grad = grad
        return grad

    def _accum_poly_grad(self, pol, g, grad, vp):
        for m, c in pol.items():
            L = len(m)
            if L == 0:
                continue
            elif L == 1:
                gv = grad[m[0]]
                if gv:
                    for k, cc in gv.items(): g[k] = (g.get(k, 0) + c * cc) % p
            elif L == 2:
                a, b = m
                if a == b:
                    fac = (2 * c * vp[a]) % p; gv = grad[a]
                    if gv and fac:
                        for k, cc in gv.items(): g[k] = (g.get(k, 0) + fac * cc) % p
                else:
                    fa = (c * vp[b]) % p; ga = grad[a]
                    if ga and fa:
                        for k, cc in ga.items(): g[k] = (g.get(k, 0) + fa * cc) % p
                    fb = (c * vp[a]) % p; gb = grad[b]
                    if gb and fb:
                        for k, cc in gb.items(): g[k] = (g.get(k, 0) + fb * cc) % p
            else:
                for idx in range(L):
                    fac = c
                    for jj, v in enumerate(m):
                        if jj != idx: fac = (fac * vp[v]) % p
                    gv = grad[m[idx]]
                    if gv and fac:
                        for k, cc in gv.items(): g[k] = (g.get(k, 0) + fac * cc) % p

    def root_grad(self, i):
        g = {}
        self._accum_poly_grad(self.root_poly[i], g, self.grad, self.valp)
        return {k: c for k, c in g.items() if c}
