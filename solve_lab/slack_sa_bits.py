#!/usr/bin/env python3
"""Two-level slack-active search: OUTER over control bits (which set the rigid
3183-slack ripple, hence R), INNER var-repair to clean up. The plain slack_sa only
repairs continuous vars from a fixed activator and cannot change the ripple; here a
fraction of moves TOGGLE a control bit and re-seed the slack-active state (full
re-eval), exploring different ripples -> the lever for hitting a self-cancelling
(R=0) bit-setting. Squares replaced by roots Q=0. Saves any full solution."""
import json, time, random, math, sys
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
    control = list(json.load(open('control_bits.json')))
    ACT = int(sys.argv[1]) if len(sys.argv) > 1 else 1858
    SEED = int(sys.argv[2]) if len(sys.argv) > 2 else 1234
    OUT = sys.argv[3] if len(sys.argv) > 3 else 'cand_bits_SOLVED.json'

    A = list(A0); A[40782] = try_sqrt(A0[40782]); A[39550] = try_sqrt(A0[39550])
    var_atoms = defaultdict(list)
    for a, poly in enumerate(A):
        for v in atom_vars(poly): var_atoms[v].append(a)
    def viol(vv): return set(a for a in range(len(A)) if atom_resid(A[a], vv) != 0)

    def seed_state(bits):
        v1 = solve(list(bestval), bits)
        if v1[12779] == 0:  # need slack gate on
            return None
        frozen = {24026: v1[18274]-v1[35186], 27116: v1[17728]-v1[1642]}
        return run(list(v1), frozen)

    curbits = [ACT]
    val = seed_state(curbits); bad = viol(val)
    best_bad = len(bad); best_val = list(val); best_bits = list(curbits)
    print(f"seed act={ACT}: {len(bad)} violated ({time.time()-t0:.0f}s)", flush=True)
    PROT = {9770, 3183, 18274, 17728, 24026, 27116, 12779, 14402}
    rng = random.Random(SEED)
    T = 3.0; it = 0; TIME = 5400
    b22 = [b for b in control if b in set([1782,1858,2795,2800,3483,5443,10652,19520,21188,21588,23634,26947,27512,29682,30104,30596,30658,30792,33251,37748,37885,38116])]
    while time.time()-t0 < TIME and best_bad > 0:
        it += 1
        # 3% of iterations: OUTER bit toggle (re-seed ripple)
        if rng.random() < 0.03:
            nb = list(curbits)
            b = control[rng.randrange(len(control))]
            if b in nb: nb.remove(b)
            else: nb.append(b)
            ns = seed_state(nb)
            if ns is None:  # ensure gate on by adding a 22-activator
                nb2 = nb + [b22[rng.randrange(len(b22))]]
                ns = seed_state(nb2)
                if ns is None: continue
                nb = nb2
            nbad = viol(ns); d = len(nbad) - len(bad)
            if d <= 0 or rng.random() < math.exp(-d / max(T, 0.05)):
                curbits = nb; val = ns; bad = nbad
        else:
            # INNER var repair
            if not bad: break
            a = rng.choice(tuple(bad)); poly = A[a]
            cand = [v for v in atom_vars(poly) if v not in PROT and v not in control]
            if not cand: bad.discard(a); continue
            v = rng.choice(cand); nv = solve_for(poly, val, v)
            if nv is None or nv == val[v]: continue
            old = val[v]
            before = sum(1 for aa in var_atoms[v] if atom_resid(A[aa], val) != 0)
            val[v] = nv
            after = sum(1 for aa in var_atoms[v] if atom_resid(A[aa], val) != 0)
            d = after - before
            if d <= 0 or rng.random() < math.exp(-d / max(T, 0.05)):
                for aa in var_atoms[v]:
                    if atom_resid(A[aa], val) != 0: bad.add(aa)
                    else: bad.discard(aa)
            else:
                val[v] = old
        if len(bad) < best_bad:
            best_bad = len(bad); best_val = list(val); best_bits = list(curbits)
            print(f"  it {it}: NEW BEST {best_bad}: {sorted(bad)} bits={sorted(curbits)[:8]} ({time.time()-t0:.0f}s)", flush=True)
            if best_bad == 0:
                allbad = viol_atoms(A0, val)
                print(f"  ORIGINAL verify: {len(allbad)}", flush=True)
                if not allbad:
                    json.dump({f"x_{i}": val[i] for i in range(NVARS)}, open(OUT, 'w'))
                    print(f"  *** SOLVED *** -> {OUT}", flush=True); return
                bad = viol(val); best_bad = len(bad)
        if it % 4000 == 0:
            T *= 0.9
            if T < 0.1: T = 3.0
        if it % 60000 == 0:
            print(f"  it {it}: bad={len(bad)} best={best_bad} T={T:.2f} |bits|={len(curbits)} ({time.time()-t0:.0f}s)", flush=True)
            val = list(best_val); curbits = list(best_bits); bad = viol(val)
    print(f"done: best {best_bad} bits={sorted(best_bits)} ({time.time()-t0:.0f}s)", flush=True)
    json.dump({"bad": sorted(viol(best_val)), "bits": sorted(best_bits), "val": {str(i): best_val[i] for i in range(NVARS)}}, open(OUT.replace('.json','_best.json'), 'w'))

if __name__ == '__main__':
    main()
