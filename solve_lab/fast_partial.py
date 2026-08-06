#!/usr/bin/env python3
"""Fast PARTIAL forward-evaluator: to compute a few target vars, only evaluate their
transitive ancestors in the confluent DAG (not all 31k vars). ~100x faster for
single-target queries -> unblocks rapid structural experiments. Returns a function
peval(setbits, targets) -> {var: value}."""
import json
from collections import deque
from confluent_eval5 import build5, make_forward

def build_partial():
    A, kind, info, seq, bestval, ncyc = build5()
    pos = {v: i for i, v in enumerate(seq)}
    # dependency map: var -> its input vars (from info)
    deps = {}
    for v in seq:
        k = kind[v]; d = set()
        if k == 'gate':
            coef, terms = info[v]
            for c, m in terms: d.update(m)
        elif k == 'load':
            bit, cbx, lt = info[v]; d.add(bit)
            for c, m in lt: d.update(m)
        elif k == 'div':
            c, u, rest = info[v]; d.add(u)
            for cc, m in rest: d.update(m)
        d.discard(v); deps[v] = d
    defined = set(seq)

    def ancestors(targets):
        seen = set(); q = deque(t for t in targets if t in defined)
        while q:
            v = q.popleft()
            if v in seen: continue
            seen.add(v)
            for u in deps.get(v, ()):
                if u in defined and u not in seen: q.append(u)
        # order by seq position
        return sorted(seen, key=lambda v: pos[v])

    def peval(setbits, targets, base=None):
        val = list(bestval) if base is None else base
        for b in setbits: val[b] = 1
        sub = ancestors(targets)
        for v in sub:
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
                    num = -rest; den = cbx*val[bit]
                    if den and num % den == 0: val[v] = num//den
            elif k == 'div':
                c, u, rest = info[v]; rs = 0
                for cc, m in rest:
                    t = cc
                    for x in m: t *= val[x]
                    rs += t
                den = c*val[u]
                if den and (-rs) % den == 0: val[v] = (-rs)//den
                elif den == 0: val[v] = 0
        return val
    return peval, ancestors, A, kind, info, bestval

if __name__ == '__main__':
    import time
    t0 = time.time()
    peval, anc, A, kind, info, bestval = build_partial()
    print(f"built ({time.time()-t0:.0f}s)", flush=True)
    # sanity: x_9770 ancestors and value
    for tv in (9770, 3183, 18274, 17728, 8821):
        a = anc([tv])
        print(f"x_{tv}: {len(a)} ancestors", flush=True)
    t1 = time.time()
    for _ in range(50):
        v = peval([1858], [9770, 3183, 18274, 17728, 8821])
    print(f"50 partial evals of 5 targets: {time.time()-t1:.2f}s  (x_9770={v[9770]})", flush=True)
