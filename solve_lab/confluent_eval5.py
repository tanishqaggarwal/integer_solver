#!/usr/bin/env python3
"""Confluent evaluator v5 — handles DIVISION-oriented wires (major fix).

v4 froze 1504 wires whose defining gate places the target inside a product
(e.g. x_18274 from 4954: x_6773 = x_8821*x_18274 -> x_18274 = x_6773/x_8821).
best solved these by division; v4's linear-only solver left them constant, which
wrongly froze x_18274 (and made the twist look unsatisfiable). Here we add:
  - load wires  : x_B = bit*(HUGE+s*x_C)   (huge-atoms)
  - div wires   : v in exactly one product (u,v): v = -rest/(c*u)
plus a proper topological order for the whole orientation. Exact Z or mod-P."""
import json, time, sys
from collections import defaultdict, deque
from propagate import load_atoms, atom_vars, NVARS

def boolean_vars(atoms):
    b = set()
    for poly in atoms:
        if len(poly) == 2:
            ms = sorted(poly.keys(), key=len)
            if len(ms[0]) == 1 and len(ms[1]) == 2 and ms[1] == (ms[0][0], ms[0][0]) and poly[ms[0]] == -poly[ms[1]]:
                b.add(ms[0][0])
    return b

def build5():
    A = load_atoms(); bset = boolean_vars(A)
    prov = json.load(open('eval_order.json'))['prov']
    best = json.load(open('best/best_partial_39019.json'))
    bestval = [0]*NVARS
    for k, x in best.items(): bestval[int(k[2:])] = x
    df = {}
    for v in range(NVARS):
        p = prov[v] if v < len(prov) else None
        if p and p[0] >= 0: df[v] = p[0]
    df[9770] = 27973; df[3183] = 27978
    # huge-atom loads
    load_atom = {}
    for a, poly in enumerate(A):
        prod2 = [(m, c) for m, c in poly.items() if len(m) == 2 and m[0] != m[1]]
        if len(prod2) != 1: continue
        (m2, c2) = prod2[0]
        big = any(abs(c) > 10**40 for m, c in poly.items() if len(m) == 1)
        for bit, xB in [(m2[0], m2[1]), (m2[1], m2[0])]:
            if bit in bset and big and xB not in load_atom:
                load_atom[xB] = (a, bit); break
    # classify each defined var
    kind = {}      # v -> ('gate'|'load'|'div'|'const')
    info = {}
    for v in df:
        a = df[v]; poly = A[a]
        if v in load_atom:
            la, bit = load_atom[v]
            cbx = next(c for m, c in A[la].items() if len(m) == 2 and v in m and bit in m)
            lterms = [(c, m) for m, c in A[la].items() if not (len(m) == 2 and v in m and bit in m)]
            kind[v] = 'load'; info[v] = (bit, cbx, lterms); continue
        mons_with_v = [(m, c) for m, c in poly.items() if v in m]
        if (v,) in poly:
            coef = poly[(v,)]; terms = [(c, m) for m, c in poly.items() if m != (v,)]
            kind[v] = 'gate'; info[v] = (coef, terms)
        elif len(mons_with_v) == 1 and len(mons_with_v[0][0]) == 2:
            (m2, c) = mons_with_v[0]
            u = m2[0] if m2[1] == v else m2[1]
            rest = [(cc, mm) for mm, cc in poly.items() if v not in mm]
            kind[v] = 'div'; info[v] = (c, u, rest)
        else:
            kind[v] = 'const'   # complex (v^2 etc.) - leave at best
    # dependencies
    deps = {}
    for v in df:
        d = set()
        if kind[v] == 'load':
            bit, cbx, lt = info[v]; d.add(bit)
            for c, m in lt: d.update(m)
        elif kind[v] == 'gate':
            coef, terms = info[v]
            for c, m in terms: d.update(m)
        elif kind[v] == 'div':
            c, u, rest = info[v]; d.add(u)
            for cc, m in rest: d.update(m)
        d.discard(v); deps[v] = d
    defined = set(v for v in df if kind[v] != 'const')
    indeg = {v: 0 for v in defined}; adj = defaultdict(list)
    for v in defined:
        for x in deps[v]:
            if x in defined: adj[x].append(v); indeg[v] += 1
    q = deque([v for v in defined if indeg[v] == 0]); topo = []
    while q:
        v = q.popleft(); topo.append(v)
        for u in adj[v]:
            indeg[u] -= 1
            if indeg[u] == 0: q.append(u)
    cyc = [v for v in defined if v not in set(topo)]
    seq = topo + cyc
    return A, kind, info, seq, bestval, len(cyc)

def make_forward(kind, info, seq, bestval, mod=None):
    if mod is None:
        def solve(val, setbits, watch=None):
            for b in setbits: val[b] = 1
            for v in seq:
                k = kind[v]
                if k == 'gate':
                    coef, terms = info[v]; rs = 0
                    for c, m in terms:
                        t = c
                        for x in m: t *= val[x]
                        rs += t
                    if coef and (-rs) % coef == 0: val[v] = (-rs)//coef
                elif k == 'load':
                    bit, cbx, lt = info[v]
                    if val[bit] == 0: val[v] = 0
                    else:
                        rest = 0
                        for c, m in lt:
                            t = c
                            for x in m: t *= (1 if x == bit else val[x])
                            rest += t
                        num = -rest; den = cbx * val[bit]
                        if den and num % den == 0: val[v] = num//den
                elif k == 'div':
                    c, u, rest = info[v]; rs = 0
                    for cc, m in rest:
                        t = cc
                        for x in m: t *= val[x]
                        rs += t
                    den = c * val[u]
                    if den and (-rs) % den == 0: val[v] = (-rs)//den
                    elif den == 0: val[v] = 0
            return val
        return solve
    P = mod
    from modp import inv
    def solveP(val, setbits, watch=None):
        for b in setbits: val[b] = 1
        for v in seq:
            k = kind[v]
            if k == 'gate':
                coef, terms = info[v]
                if coef % P == 0: continue
                rs = 0
                for c, m in terms:
                    t = c % P
                    for x in m: t = (t*val[x]) % P
                    rs = (rs+t) % P
                val[v] = (-rs * inv(coef)) % P
            elif k == 'load':
                bit, cbx, lt = info[v]
                if val[bit] == 0: val[v] = 0
                else:
                    rest = 0
                    for c, m in lt:
                        t = c % P
                        for x in m: t = (t*(1 if x == bit else val[x])) % P
                        rest = (rest+t) % P
                    val[v] = (-rest * inv((cbx*val[bit]) % P)) % P
            elif k == 'div':
                c, u, rest = info[v]
                if val[u] == 0: val[v] = 0; continue
                rs = 0
                for cc, m in rest:
                    t = cc % P
                    for x in m: t = (t*val[x]) % P
                    rs = (rs+t) % P
                val[v] = (-rs * inv((c*val[u]) % P)) % P
        return val
    return solveP

def main():
    t0 = time.time()
    A, kind, info, seq, bestval, ncyc = build5()
    nk = defaultdict(int)
    for v in kind: nk[kind[v]] += 1
    print(f"kinds: {dict(nk)}, seq {len(seq)}, cyclic {ncyc} ({time.time()-t0:.0f}s)", flush=True)
    solve = make_forward(kind, info, seq, bestval)
    val = solve(list(bestval), [])
    viol = []
    for a, poly in enumerate(A):
        s = 0
        for m, c in poly.items():
            t = c
            for x in m: t *= val[x]
            s += t
        if s: viol.append(a)
    print(f"forward_Z([]) violated: {len(viol)} {sorted(viol)[:8]}", flush=True)
    print("VALID (==39019)" if sorted(viol) == [1817, 30378, 40782, 44271] else f"CHANGED: {sorted(viol)}", flush=True)

if __name__ == '__main__':
    main()
