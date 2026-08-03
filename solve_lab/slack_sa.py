#!/usr/bin/env python3
"""Simulated-annealing repair in SLACK-ACTIVE space, with the two degree-4 square
atoms (a40782, a39550) REPLACED by their degree-2 roots Q=0 (extracted via
try_sqrt). This makes them solvable in the min-conflicts step (previously they were
unfixable, causing the plateau at 9). SA accepts worse moves via temperature to
escape the frustrated core. Saves any full solution (verified against ORIGINAL
atoms)."""
import json, time, random, math
from collections import defaultdict
from confluent_eval5 import build5, make_forward
from slack_active import make_slack_solver, viol_atoms
from check_square import try_sqrt
from propagate import atom_vars, NVARS

def atom_resid(poly, val):
    s = 0
    for m, c in poly.items():
        t = c
        for x in m: t *= val[x]
        s += t
    return s

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
        s = math.isqrt(r)
        return s if s*s == r else None
    # general quadratic c2 v^2 + c1 v + c0 = 0
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
    control = list(json.load(open('control_bits.json')))

    # build A' with squares replaced by roots
    SQ = {40782: try_sqrt(A0[40782]), 39550: try_sqrt(A0[39550])}
    A = list(A0)
    for a, Q in SQ.items():
        A[a] = Q
    var_atoms = defaultdict(list)
    for a, poly in enumerate(A):
        for v in atom_vars(poly): var_atoms[v].append(a)

    def viol(vv):
        return set(a for a in range(len(A)) if atom_resid(A[a], vv) != 0)

    # seed from bit 1858 slack-active
    v1 = solve(list(bestval), [1858])
    frozen = {24026: v1[18274]-v1[35186], 27116: v1[17728]-v1[1642]}
    val = run(list(v1), frozen)
    bad = viol(val)
    print(f"seed: {len(bad)} violated (A' with roots) ({time.time()-t0:.0f}s)", flush=True)

    PROT = {9770, 3183, 18274, 17728, 24026, 27116, 12779, 14402}
    rng = random.Random(1234)
    best_bad = len(bad); best_val = list(val)
    T = 3.0; it = 0
    TIME = 3300
    while time.time()-t0 < TIME and best_bad > 0:
        it += 1
        if not bad:
            break
        a = rng.choice(tuple(bad))
        poly = A[a]
        cand = [v for v in atom_vars(poly) if v not in PROT]
        if not cand:
            bad.discard(a); continue
        v = rng.choice(cand)
        nv = solve_for(poly, val, v)
        if nv is None or nv == val[v]:
            continue
        old = val[v]
        before = sum(1 for aa in var_atoms[v] if atom_resid(A[aa], val) != 0)
        val[v] = nv
        after = sum(1 for aa in var_atoms[v] if atom_resid(A[aa], val) != 0)
        delta = after - before
        # SA acceptance
        if delta <= 0 or rng.random() < math.exp(-delta / max(T, 0.05)):
            for aa in var_atoms[v]:
                if atom_resid(A[aa], val) != 0: bad.add(aa)
                else: bad.discard(aa)
        else:
            val[v] = old
        if len(bad) < best_bad:
            best_bad = len(bad); best_val = list(val)
            print(f"  it {it}: NEW BEST {best_bad}: {sorted(bad)} (T={T:.2f}, {time.time()-t0:.0f}s)", flush=True)
            if best_bad == 0:
                allbad = viol_atoms(A0, val)  # verify ORIGINAL atoms
                print(f"  ORIGINAL-atom verify: {len(allbad)} violated", flush=True)
                if not allbad:
                    json.dump({f"x_{i}": val[i] for i in range(NVARS)}, open('cand_SA_SOLVED.json','w'))
                    print("  *** SOLVED ***", flush=True); return
                bad = viol(val); best_bad = len(bad)
        # cooling / reheating
        if it % 4000 == 0:
            T *= 0.85
            if T < 0.1: T = 3.0  # reheat
        if it % 40000 == 0:
            print(f"  it {it}: bad={len(bad)} best={best_bad} T={T:.2f} ({time.time()-t0:.0f}s)", flush=True)
            # restart from best occasionally
            val = list(best_val); bad = viol(val)
    print(f"SA done: best {best_bad} ({time.time()-t0:.0f}s)", flush=True)
    json.dump({"bad": sorted(viol(best_val)), "val": {str(i): best_val[i] for i in range(NVARS)}}, open('slack_sa_best.json','w'))

if __name__ == '__main__':
    main()
