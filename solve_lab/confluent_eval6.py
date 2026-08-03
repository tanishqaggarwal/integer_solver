#!/usr/bin/env python3
"""v6 — product-FORWARD orientation (cascade fix attempt).

The v5 cascade is orientation: best defines products' outputs via other gates, so
when a product input changes the output doesn't track (product atom violates). v6
forces every product/square/load to DEFINE its output var (forward). Sums fill in
the rest; combos are checks. If acyclic and it reproduces 39,019 at all-0, then a
bit flip's residue load propagates through products correctly (no product cascade).
"""
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

def build6():
    A = load_atoms(); bset = boolean_vars(A)
    prov = json.load(open('eval_order.json'))['prov']
    best = json.load(open('best/best_partial_39019.json'))
    bestval = [0]*NVARS
    for k, x in best.items(): bestval[int(k[2:])] = x
    NV_PRIM = 4
    prim = [len(atom_vars(A[a])) <= NV_PRIM for a in range(len(A))]
    # classify primitive atoms
    load_atom = {}
    for a, poly in enumerate(A):
        prod2 = [(m, c) for m, c in poly.items() if len(m) == 2 and m[0] != m[1]]
        if len(prod2) != 1: continue
        (m2, c2) = prod2[0]
        big = any(abs(c) > 10**40 for m, c in poly.items() if len(m) == 1)
        for bit, xB in [(m2[0], m2[1]), (m2[1], m2[0])]:
            if bit in bset and big and xB not in load_atom:
                load_atom[xB] = (a, bit); break
    # linear-output candidates
    def lin_outs(a):
        bad = set(); lin = set()
        for m in A[a]:
            if len(m) == 1: lin.add(m[0])
            else: bad.update(m)
        return lin - bad
    has_prod = [any(len(m) >= 2 for m in A[a]) for a in range(len(A))]
    kind = {}; info = {}
    match = {}; used = [False]*len(A)   # var->atom, atom used
    # 1) loads (forced)
    for xB, (a, bit) in load_atom.items():
        if xB not in match and not used[a]:
            match[xB] = ('load', a, bit); used[a] = True
    # 2) products/squares forward: output = the linear var
    for a in range(len(A)):
        if not prim[a] or used[a] or not has_prod[a]: continue
        outs = lin_outs(a)
        for v in outs:
            if v not in match:
                match[v] = ('prodfwd', a, None); used[a] = True; break
    # 3) sums/affines (linear primitives) -> remaining var
    def augment(a, seenv):
        for v in lin_outs(a):
            if v in seenv: continue
            seenv.add(v)
            cur = match.get(v)
            if cur is None:
                match[v] = ('gate', a, None); used[a] = True; return True
            if cur[0] in ('load', 'prodfwd'): continue  # don't displace forced
            if augment(cur[1], seenv):
                match[v] = ('gate', a, None); used[a] = True; return True
        return False
    sys.setrecursionlimit(1_000_000)
    for a in range(len(A)):
        if not prim[a] or used[a] or has_prod[a]: continue
        augment(a, set())
    # 4) division fallback for still-unmatched vars best defined via a product
    for v in range(NVARS):
        if v in match: continue
        p = prov[v] if v < len(prov) else None
        if not p or p[0] < 0: continue
        a = p[0]; poly = A[a]
        if used[a]: continue
        mons_with_v = [(m, c) for m, c in poly.items() if v in m]
        if (v,) not in poly and len(mons_with_v) == 1 and len(mons_with_v[0][0]) == 2 and len(atom_vars(poly)) <= NV_PRIM:
            match[v] = ('div', a, None); used[a] = True
    # compile per-var evaluators
    for v, (k, a, extra) in match.items():
        poly = A[a]
        if k == 'load':
            bit = extra
            cbx = next(c for m, c in poly.items() if len(m) == 2 and v in m and bit in m)
            lterms = [(c, m) for m, c in poly.items() if not (len(m) == 2 and v in m and bit in m)]
            kind[v] = 'load'; info[v] = (bit, cbx, lterms)
        elif k in ('prodfwd', 'gate'):
            coef = poly.get((v,), 0); terms = [(c, m) for m, c in poly.items() if m != (v,)]
            kind[v] = 'gate'; info[v] = (coef, terms)
        elif k == 'div':
            mons_with_v = [(m, c) for m, c in poly.items() if v in m]
            (m2, c) = mons_with_v[0]; u = m2[0] if m2[1] == v else m2[1]
            rest = [(cc, mm) for mm, cc in poly.items() if v not in mm]
            kind[v] = 'div'; info[v] = (c, u, rest)
    # deps + topo
    deps = {}
    for v in kind:
        d = set()
        if kind[v] == 'load':
            bit, cbx, lt = info[v]; d.add(bit)
            for c, m in lt: d.update(m)
        elif kind[v] == 'gate':
            for c, m in info[v][1]: d.update(m)
        elif kind[v] == 'div':
            c, u, rest = info[v]; d.add(u)
            for cc, m in rest: d.update(m)
        d.discard(v); deps[v] = d
    defined = set(kind)
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

def main():
    t0 = time.time()
    A, kind, info, seq, bestval, ncyc = build6()
    nk = defaultdict(int)
    for v in kind: nk[kind[v]] += 1
    print(f"v6 kinds {dict(nk)}, defined {len(kind)}, cyclic {ncyc} ({time.time()-t0:.0f}s)", flush=True)
    val = list(bestval)
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
                den = cbx*val[bit]
                if den and (-rest) % den == 0: val[v] = (-rest)//den
        elif k == 'div':
            c, u, rest = info[v]
            if val[u] == 0: val[v] = 0; continue
            rs = 0
            for cc, m in rest:
                t = cc
                for x in m: t *= val[x]
                rs += t
            den = c*val[u]
            if den and (-rs) % den == 0: val[v] = (-rs)//den
    viol = []
    for a, poly in enumerate(A):
        s = 0
        for m, c in poly.items():
            t = c
            for x in m: t *= val[x]
            s += t
        if s: viol.append(a)
    print(f"v6 forward([]) violated: {len(viol)} {sorted(viol)[:10]}", flush=True)
    print("VALID==39019" if sorted(viol) == [1817,30378,40782,44271] else f"different: {sorted(viol)[:15]}", flush=True)

if __name__ == '__main__':
    main()
