#!/usr/bin/env python3
"""Follow x_15690 (numerator core of x_18274) and x_14494 (of x_17728) down their
definition tree. At each wire test mod-P linearity in the 233 bits. Find the
shallowest NONLINEAR wire and characterize it (is it a product of two linear forms?
a square? a division?). If nonlinearity is one shallow product L1*L2, we can attack
via factoring; if pervasive, it's hopeless."""
import json
from collections import deque
from confluent_eval5 import build5, make_forward
P = (1 << 61) - 1
BITS22 = set([1782,1858,2795,2800,3483,5443,10652,19520,21188,21588,23634,26947,27512,29682,30104,30596,30658,30792,33251,37748,37885,38116])

def main():
    A, kind, info, seq0, bestval, ncyc = build5()
    order = json.load(open('eval_order.json'))['order']
    defset = set(v for v in kind if kind[v] != 'const')
    seq = [v for v in order if v in defset and v not in (9770, 3183)]
    seq += [v for v in (9770, 3183) if v in defset]
    seq += [v for v in defset if v not in set(order) and v not in (9770, 3183)]
    solveP = make_forward(kind, info, seq, bestval, mod=P)
    bm = [x % P for x in bestval]
    control = json.load(open('control_bits.json'))
    bits233 = [b for b in control if b not in BITS22]

    b0 = solveP(list(bm), [])
    single = {b: solveP(list(bm), [b]) for b in bits233}
    def is_linear(v):
        for sizes, seedmul in [([3,10,40,120], 71), ([5,25,90], 53)]:
            for si, sz in enumerate(sizes):
                S = sorted(set(bits233[(i*seedmul+si*131+3) % len(bits233)] for i in range(sz)))
                val = solveP(list(bm), S)
                pred = (b0[v] + sum((single[b][v]-b0[v]) for b in S)) % P
                if val[v] != pred: return False
        return True

    def deps_of(v):
        k = kind.get(v); d = []
        if k == 'gate':
            for c, m in info[v][1]:
                for x in m: d.append((x, len(m)))
        elif k == 'div':
            c, u, rest = info[v]; d.append((u, 1))
            for cc, m in rest:
                for x in m: d.append((x, len(m)))
        elif k == 'load':
            bit, cbx, lt = info[v]
            for c, m in lt:
                for x in m: d.append((x, len(m)))
        return d

    for root in (15690, 14494):
        print(f"\n===== tracing x_{root} (linear={is_linear(root)}) =====")
        # BFS to find shallowest nonlinear wires and structure
        seen = set(); q = deque([(root, 0)])
        shallow_nl = []
        while q and len(shallow_nl) < 12:
            v, depth = q.popleft()
            if v in seen or depth > 25: continue
            seen.add(v)
            k = kind.get(v)
            lin = is_linear(v) if k else True
            if not lin and k:
                # characterize: does it have a product term (two non-bit vars)?
                if k == 'gate':
                    coef, terms = info[v]
                    prods = [(c, m) for c, m in terms if len(m) >= 2]
                    tag = f"gate {len(terms)}terms prods={len(prods)}"
                    if prods:
                        m = prods[0][1]
                        tag += f" e.g. x_{m[0]}(lin={is_linear(m[0])})*x_{m[1] if len(m)>1 else m[0]}(lin={is_linear(m[1] if len(m)>1 else m[0])})"
                elif k == 'div':
                    c, u, rest = info[v]; tag = f"div u=x_{u}(lin={is_linear(u)})"
                elif k == 'load':
                    tag = "load"
                else: tag = k
                shallow_nl.append((v, depth, tag))
            for x, mdeg in deps_of(v):
                if x not in seen: q.append((x, depth+1))
        for v, depth, tag in shallow_nl:
            print(f"  NONLIN x_{v} @depth {depth}: {tag}", flush=True)
    print("done", flush=True)

if __name__ == '__main__':
    main()
