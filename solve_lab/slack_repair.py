#!/usr/bin/env python3
"""Min-conflicts integer repair seeded from the SLACK-ACTIVE state. The twist
holds by construction; ~18 secondary atoms are broken because activating x_12779=1
perturbs coupled wires. Each broken atom is simple; repair by solving it for one of
its variables (linear -> divide; product -> divide), WalkSAT-style: pick a random
broken atom, try each 'solve-for' variable, keep the move that most reduces total
violations. Runs long; saves any full solution."""
import json, time, random
from collections import defaultdict
from confluent_eval5 import build5, make_forward
from slack_active import make_slack_solver, viol_atoms
from propagate import atom_vars, NVARS

def atom_resid(poly, val):
    s = 0
    for m, c in poly.items():
        t = c
        for x in m: t *= val[x]
        s += t
    return s

def solve_for(poly, val, v):
    """Solve poly==0 for variable v given others in val. Returns int or None."""
    # collect coefficient of v across monomials: poly = A*v + B  (if v linear),
    # or product terms with v. Group by power of v.
    c0 = 0; c1 = 0; c2 = 0  # constant, coef*v, coef*v^2
    for m, c in poly.items():
        k = m.count(v)
        t = c
        for x in m:
            if x != v: t *= val[x]
        if k == 0: c0 += t
        elif k == 1: c1 += t
        else: c2 += t  # v^2 * (rest)
    if c2 == 0:
        if c1 == 0: return None
        if (-c0) % c1 != 0: return None
        return (-c0)//c1
    # c2*v^2 + c1*v + c0 = 0
    if c1 == 0:
        if c2 == 0: return None
        q = -c0
        if q % c2 != 0: return None
        r = q//c2
        if r < 0: return None
        import math
        s = math.isqrt(r)
        if s*s == r: return s
        return None
    return None

def main():
    t0 = time.time()
    A, kind, info, seq0, bestval, ncyc = build5()
    order = json.load(open('eval_order.json'))['order']
    defset = set(v for v in kind if kind[v] != 'const')
    seq = [v for v in order if v in defset and v not in (9770, 3183)]
    seq += [v for v in (9770, 3183) if v in defset]
    seq += [v for v in defset if v not in set(order) and v not in (9770, 3183)]
    solve = make_forward(kind, info, seq, bestval)
    run, seq2 = make_slack_solver(kind, info, seq, bestval)
    control = set(json.load(open('control_bits.json')))
    pins = set(json.load(open('pins.json'))) if False else set()

    # var -> atoms index
    var_atoms = defaultdict(list)
    for a, poly in enumerate(A):
        for v in atom_vars(poly):
            var_atoms[v].append(a)

    # seed slack-active state from best activator bit 1858
    v1 = solve(list(bestval), [1858])
    frozen = {24026: v1[18274]-v1[35186], 27116: v1[17728]-v1[1642]}
    val = run(list(v1), frozen)
    bad = set(viol_atoms(A, val))
    print(f"seed slack-active: {len(bad)} violated ({time.time()-t0:.0f}s)", flush=True)

    # protected vars: never touch (twist targets + slack inputs + pins)
    PROT = {9770, 3183, 18274, 17728, 24026, 27116, 12779, 14402}
    rng = random.Random(7)
    best_bad = len(bad)
    it = 0
    # precompute residual cache
    def total_bad(vv):
        return set(viol_atoms(A, vv))

    while time.time()-t0 < 3000 and best_bad > 0:
        it += 1
        if not bad: break
        a = rng.choice(list(bad))
        poly = A[a]
        cand = [v for v in atom_vars(poly) if v not in PROT]
        if not cand:
            bad.discard(a); continue
        best_move = None; best_delta = 1  # require strict improvement (or equal w/ prob)
        rng.shuffle(cand)
        for v in cand[:12]:
            old = val[v]
            nv = solve_for(poly, val, v)
            if nv is None or nv == old: continue
            # evaluate local delta: atoms touching v
            before = sum(1 for aa in var_atoms[v] if atom_resid(A[aa], val) != 0)
            val[v] = nv
            after = sum(1 for aa in var_atoms[v] if atom_resid(A[aa], val) != 0)
            delta = after - before
            if delta < best_delta or (delta == best_delta and rng.random() < 0.3):
                best_delta = delta; best_move = (v, nv, old)
            val[v] = old
        if best_move is not None:
            v, nv, old = best_move
            val[v] = nv
            # update bad set locally
            for aa in var_atoms[v]:
                if atom_resid(A[aa], val) != 0: bad.add(aa)
                else: bad.discard(aa)
        else:
            # random walk: pick any solvable var, jump
            v = rng.choice(cand)
            nv = solve_for(poly, val, v)
            if nv is not None:
                val[v] = nv
                for aa in var_atoms[v]:
                    if atom_resid(A[aa], val) != 0: bad.add(aa)
                    else: bad.discard(aa)
        if len(bad) < best_bad:
            best_bad = len(bad)
            print(f"  it {it}: NEW BEST {best_bad} violated: {sorted(bad)} ({time.time()-t0:.0f}s)", flush=True)
            json.dump({"bad": sorted(bad), "val": {str(i): val[i] for i in range(NVARS)}}, open('slack_repair_best.json','w'))
            if best_bad == 0:
                # full verify
                allbad = total_bad(val)
                print(f"  full verify: {len(allbad)} violated", flush=True)
                if not allbad:
                    json.dump({f"x_{i}": val[i] for i in range(NVARS)}, open('cand_repair_SOLVED.json','w'))
                    print("  *** SOLVED (slack-active repair) ***", flush=True); return
                else:
                    bad = allbad; best_bad = len(bad)
        if it % 2000 == 0:
            print(f"  it {it}: bad={len(bad)} best={best_bad} ({time.time()-t0:.0f}s)", flush=True)
    print(f"repair done: best {best_bad} violated after {it} iters ({time.time()-t0:.0f}s)", flush=True)

if __name__ == '__main__':
    main()
