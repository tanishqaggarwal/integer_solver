#!/usr/bin/env python3
"""Re-orient x_18274 and x_17728 away from their DIV wires (a4954: x_8821*x_18274=
x_6773; a13204: x_8821*x_17728=x_17233) -- which quantize them to g2*Z/h2*Z and
force degeneracy -- onto their LINEAR GATE atoms:
  x_18274 = x_31434 - x_6283      (a11398)
  x_17728 = x_27912 - x_28035     (a11388)
The div atoms a4954/a13204 become CHECKS (x_6773 = x_8821*x_18274, etc.). This lets
x_18274/x_17728 escape quantization. Rebuild the topo order and test:
 (1) forward_Z([]) validity, (2) does x_18274 now escape g2*Z, (3) can x_18274 reach
 a x_9770 (22-side) value / x_17728 reach a x_3183 value."""
import json, time
from collections import defaultdict, deque
from math import gcd
from confluent_eval5 import build5, make_forward
from propagate import NVARS

BITS22=set([1782,1858,2795,2800,3483,5443,10652,19520,21188,21588,23634,26947,27512,29682,30104,30596,30658,30792,33251,37748,37885,38116])

def rebuild_with_overrides(A, kind, info, overrides):
    """overrides: var -> defining atom index (must be a linear 'gate' for that var).
    Recompute kind/info for the overridden vars and their displaced div-atoms."""
    kind = dict(kind); info = dict(info)
    for v, a in overrides.items():
        poly = A[a]
        assert (v,) in poly, f"atom {a} not linear in x_{v}"
        coef = poly[(v,)]
        terms = [(c, m) for m, c in poly.items() if m != (v,)]
        kind[v] = 'gate'; info[v] = (coef, terms)
    return kind, info

def make_seq(kind, info):
    defined = set(v for v in kind if kind[v] != 'const')
    deps = {}
    for v in defined:
        d = set()
        k = kind[v]
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
    return topo + cyc, len(cyc)

def main():
    t0 = time.time()
    A, kind, info, seq0, bestval, ncyc = build5()
    # override x_18274 -> a11398 (x_31434 - x_6283 + x_18274 == 0 => x_18274 = x_6283 - x_31434? check sign)
    # a11398: 1*x6283 + -1*x18274 + 1*x31434  => -x_18274 = -x6283 - x31434 => x_18274 = x6283 + x31434?
    # solve: -1*x18274 + (x6283 + x31434) = 0 => x18274 = x6283 + x31434
    # a11388: 1*x17728 + -1*x27912 + 1*x28035 => x17728 = x27912 - x28035
    over = {18274: 11398, 17728: 11388}
    kind2, info2 = rebuild_with_overrides(A, kind, info, over)
    seq2, ncyc2 = make_seq(kind2, info2)
    print(f"re-oriented seq len {len(seq2)}, cyclic {ncyc2} (was {ncyc}) ({time.time()-t0:.0f}s)", flush=True)
    solve2 = make_forward(kind2, info2, seq2, bestval)
    val = solve2(list(bestval), [])
    # validity
    viol = []
    for a, poly in enumerate(A):
        s = 0
        for m, c in poly.items():
            t = c
            for x in m: t *= val[x]
            s += t
        if s: viol.append(a)
    print(f"forward_Z([]) violated: {len(viol)}: {sorted(viol)[:12]}", flush=True)

    control = json.load(open('control_bits.json'))
    b233 = [b for b in control if b not in BITS22]
    b22 = [b for b in control if b in BITS22]
    base = solve2(list(bestval), [])
    g2 = base[18274]; g = base[9770]; h2 = base[17728]; h = base[3183]
    print(f"base x_18274={g2}, x_9770={g}, x_17728={h2}, x_3183={h}", flush=True)
    # does x_18274 escape g2*Z now?
    st = 11
    def rnd():
        nonlocal st; st=(st*6364136223846793005+1442695040888963407)&((1<<64)-1); return st>>33
    esc = 0; imgs = set()
    g2b = g2 if g2 else 1
    for _ in range(300):
        k = 1+rnd()%14
        S = sorted(set(b233[rnd()%len(b233)] for _ in range(k)))
        v = solve2(list(bestval), S)
        imgs.add(v[18274])
        if g2 and v[18274] % g2 != 0: esc += 1
    print(f"x_18274 escaped g2*Z: {esc}/300; distinct x_18274 values: {len(imgs)}", flush=True)
    print(f"done ({time.time()-t0:.0f}s)", flush=True)

if __name__ == '__main__':
    main()
