#!/usr/bin/env python3
"""CORE: trace the ESCAPE CASCADE. For x_24026 to be nonzero (slack active), a chain
of gates must be nonzero. Follow the 'must be nonzero' requirement backward through
the confluent defining atoms until we hit ROOTS: control bits, free vars (df=None),
or 0-constants. If a root is a controllable bit/free var, the cascade IS activatable
(and the whole 'trapdoor' is an artifact of forward-eval's zero-fixpoint). If it's
circular or grounded at a structural 0, it's genuinely blocked."""
import json
from collections import deque
from confluent_eval5 import build5
from propagate import load_atoms, atom_vars

def main():
    A, kind, info, seq0, bestval, ncyc = build5()
    prov = json.load(open('eval_order.json'))['prov']
    control = set(json.load(open('control_bits.json')))
    def dfa(v):
        p = prov[v] if v < len(prov) else None
        return p[0] if p and p[0] >= 0 else None

    def describe(v):
        d = dfa(v); k = kind.get(v, 'const')
        if v in control: return f"x_{v}[BIT]"
        if d is None: return f"x_{v}[FREE df=None]"
        return f"x_{v}[{k} a{d}]"

    # For a var that must be NONZERO, return the set of vars that must be nonzero
    # (for products: all factors; for linear diff a-b: 'a!=b', we flag as branch)
    def requires_nonzero(v):
        d = dfa(v)
        if d is None or v in control:
            return None  # root
        poly = A[d]
        # find the monomial(s) defining v: v is 'output'. Solve poly=0 for v.
        # gate: v linear -> v = -rest/coef ; nonzero iff rest nonzero
        # div: v = -rest/(c*u) ; nonzero iff rest nonzero
        # collect the 'rest' (terms without v)
        rest = {m: c for m, c in poly.items() if v not in m}
        rv = set()
        for m in rest:
            rv.update(m)
        return ('rest_nonzero', rest, rv)

    # BFS the cascade from x_24026
    print("=== escape cascade from x_24026 (must be nonzero) ===")
    seen = set(); q = deque([24026]); depth = {24026: 0}
    roots = []; order = []
    while q:
        v = q.popleft()
        if v in seen: continue
        seen.add(v); order.append(v)
        r = requires_nonzero(v)
        if r is None:
            roots.append(v); continue
        _, rest, rv = r
        # show the defining relation compactly
        terms = ' + '.join(('*'.join('x'+str(x) for x in m) if m else '1')+('*'+str(c) if abs(c)<10**8 else '*H') for m,c in sorted(rest.items(),key=lambda kv:-len(kv[0]))[:5])
        print(f"  {describe(v)} (depth {depth[v]}) <- needs nonzero: {terms[:90]}")
        for u in rv:
            if u not in depth: depth[u] = depth[v]+1
            if u not in seen and depth[u] <= 8: q.append(u)
    print(f"\ncascade touched {len(order)} vars")
    bit_roots = [v for v in roots if v in control]
    free_roots = [v for v in roots if dfa(v) is None and v not in control]
    print(f"ROOTS: {len(bit_roots)} control-bits, {len(free_roots)} free-vars, others {len(roots)-len(bit_roots)-len(free_roots)}")
    print(f"  bit roots: {sorted(bit_roots)[:15]}")
    print(f"  free roots: {sorted(free_roots)[:15]}")
    # key intermediate gates on the chain
    for v in [24026, 38215, 37917, 30077, 27116, 29437, 7815, 31807]:
        print(f"  {describe(v)}: value@best={bestval[v] if v<len(bestval) else '?'}")

if __name__ == '__main__':
    main()
