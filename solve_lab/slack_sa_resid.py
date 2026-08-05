#!/usr/bin/env python3
"""Slack-active SA with a RESIDUAL-BIT-LENGTH objective (not atom count). The
verifier squares have HUGE residuals (~10^250); atom-count treats them all-or-
nothing, giving the search no gradient. Here the local-move objective sums
bit_length(|residual|) over affected atoms, rewarding shrinking a huge residual
even before it hits 0 -> gradient toward R=0. Squares replaced by roots Q=0."""
import json, time, random, math, sys
from collections import defaultdict
from confluent_eval5 import build5, make_forward
from slack_active import make_slack_solver, viol_atoms
from check_square import try_sqrt
from propagate import atom_vars, NVARS

def resid(poly, val):
    s = 0
    for m, c in poly.items():
        t = c
        for x in m: t *= val[x]
        s += t
    return s

def cost(poly, val):
    r = resid(poly, val)
    return r.bit_length() if r else 0  # 0 if satisfied, else bit-length

def solve_for(poly, val, v):
    c0 = c1 = c2 = 0
    for m, c in poly.items():
        k = m.count(v); t = c
        for x in m:
            if x != v: t *= val[x]
        if k == 0: c0 += t
        elif k == 1: c1 += t
        else: c2 += t
    if c2 == 0:
        if c1 == 0 or (-c0) % c1: return None
        return (-c0)//c1
    if c1 == 0:
        if (-c0) % c2: return None
        r = (-c0)//c2
        if r < 0: return None
        s = math.isqrt(r); return s if s*s == r else None
    disc = c1*c1 - 4*c2*c0
    if disc < 0: return None
    s = math.isqrt(disc)
    if s*s != disc: return None
    for num in (-c1+s, -c1-s):
        if num % (2*c2) == 0: return num//(2*c2)
    return None

def main():
    t0 = time.time()
    A0, kind, info, seq0, bestval, ncyc = build5()
    order = json.load(open('eval_order.json'))['order']
    defset = set(v for v in kind if kind[v] != 'const')
    seq = [v for v in order if v in defset and v not in (9770, 3183)]
    seq += [v for v in (9770, 3183) if v in defset]
    seq += [v for v in defset if v not in set(order) and v not in (9770, 3183)]
    solve = make_forward(kind, info, seq, bestval)
    run, seq2 = make_slack_solver(kind, info, seq, bestval)
    ACT = int(sys.argv[1]) if len(sys.argv) > 1 else 1858
    SEED = int(sys.argv[2]) if len(sys.argv) > 2 else 55
    OUT = sys.argv[3] if len(sys.argv) > 3 else 'cand_resid_SOLVED.json'

    A = list(A0); A[40782] = try_sqrt(A0[40782]); A[39550] = try_sqrt(A0[39550])
    var_atoms = defaultdict(list)
    for a, poly in enumerate(A):
        for v in atom_vars(poly): var_atoms[v].append(a)
    def nbad(vv): return sum(1 for a in range(len(A)) if resid(A[a], vv))

    v1 = solve(list(bestval), [ACT])
    frozen = {24026: v1[18274]-v1[35186], 27116: v1[17728]-v1[1642]}
    val = run(list(v1), frozen)
    bad = set(a for a in range(len(A)) if resid(A[a], val))
    print(f"seed act={ACT}: {len(bad)} violated ({time.time()-t0:.0f}s)", flush=True)
    PROT = {9770, 3183, 18274, 17728, 24026, 27116, 12779, 14402}
    rng = random.Random(SEED)
    best_bad = len(bad); T = 40.0; it = 0; TIME = 5400
    while time.time()-t0 < TIME and best_bad > 0:
        it += 1
        if not bad: break
        a = rng.choice(tuple(bad)); poly = A[a]
        cand = [v for v in atom_vars(poly) if v not in PROT]
        if not cand: bad.discard(a); continue
        v = rng.choice(cand); nv = solve_for(poly, val, v)
        if nv is None or nv == val[v]: continue
        old = val[v]
        before = sum(cost(A[aa], val) for aa in var_atoms[v])
        val[v] = nv
        after = sum(cost(A[aa], val) for aa in var_atoms[v])
        d = after - before   # change in total bit-length cost
        if d <= 0 or rng.random() < math.exp(-d / max(T, 1.0)):
            for aa in var_atoms[v]:
                if resid(A[aa], val): bad.add(aa)
                else: bad.discard(aa)
        else:
            val[v] = old
        if len(bad) < best_bad:
            best_bad = len(bad)
            print(f"  it {it}: NEW BEST {best_bad}: {sorted(bad)} (T={T:.1f}, {time.time()-t0:.0f}s)", flush=True)
            if best_bad == 0:
                allbad = viol_atoms(A0, val)
                print(f"  ORIGINAL verify: {len(allbad)}", flush=True)
                if not allbad:
                    json.dump({f"x_{i}": val[i] for i in range(NVARS)}, open(OUT, 'w'))
                    print(f"  *** SOLVED *** -> {OUT}", flush=True); return
                bad = set(a for a in range(len(A)) if resid(A[a], val)); best_bad = len(bad)
        if it % 4000 == 0:
            T *= 0.9
            if T < 2: T = 40.0
        if it % 50000 == 0:
            print(f"  it {it}: bad={len(bad)} best={best_bad} T={T:.1f} ({time.time()-t0:.0f}s)", flush=True)
    print(f"done: best {best_bad} ({time.time()-t0:.0f}s)", flush=True)

if __name__ == '__main__':
    main()
