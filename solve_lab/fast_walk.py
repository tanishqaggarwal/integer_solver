#!/usr/bin/env python3
"""Fast mod-P DAG evaluator + WalkSAT search over the core bits.

The topological DAG evaluator computes every wire once per pass (a few passes to
settle broken cycles), all mod P so no big-int blow-up. Its violated-atom count
is 0 exactly at a true solution (the witness satisfies every definitional atom),
so minimizing it is a valid search — and it is ~20x faster than full propagation,
enabling a real WalkSAT campaign. Verifies any 0 candidate over Z."""
import json, time, sys, random
from collections import defaultdict
from fast_eval import build, NVARS
from propagate import load_atoms
from repair import boolean_vars

P = (1 << 61) - 1
_IC = {}
def inv(a):
    a %= P
    r = _IC.get(a)
    if r is None:
        r = pow(a, P - 2, P); _IC[a] = r
    return r
SEED = int(sys.argv[2]) if len(sys.argv) > 2 else 1
random.seed(SEED)

def main():
    budget = int(sys.argv[1]) if len(sys.argv) > 1 else 6000
    t0 = time.time()
    atoms, bset, define, order, base_forced, ref_val = build()
    mainv = set(json.load(open('main_comp.json'))['main_vars'])
    control = json.load(open('control_bits.json'))
    # precompute main-component atoms as mod-P monomial lists
    matoms = []
    for poly in atoms:
        vs = set()
        for m in poly: vs.update(m)
        if vs & mainv:
            matoms.append([(c % P, tuple(m)) for m, c in poly.items()])
    # precompute each wire's defining atom as monomial list (for the target solve)
    defpoly = {}
    for v in range(NVARS):
        k = define[v]
        if k[0] in ('atom', 'huge'):
            defpoly[v] = [(c % P, tuple(m)) for m, c in atoms[k[1]].items()]
    print(f"DAG+precompute in {time.time()-t0:.1f}s; main atoms {len(matoms)}", flush=True)

    def evaluate(oneset, passes=3):
        val = [0] * NVARS
        for v in order:
            if define[v][0] == 'input':
                val[v] = 1 if (v in bset and v in oneset) else (base_forced.get(v, 0) % P if v in base_forced else 0)
        # base_forced overrides (pins)
        for v, x in base_forced.items():
            val[v] = x % P
        for _ in range(passes):
            for v in order:
                dp = defpoly.get(v)
                if dp is None: continue
                c0 = c1 = 0
                for coef, m in dp:
                    if v in m:
                        # count multiplicity of v
                        t = coef; deg = 0
                        for x in m:
                            if x == v: deg += 1
                            else: t = (t * val[x]) % P
                        if deg == 1: c1 = (c1 + t) % P
                        else: c1 = c1  # skip v^2 target (rare)
                    else:
                        t = coef
                        for x in m: t = (t * val[x]) % P
                        c0 = (c0 + t) % P
                if c1:
                    val[v] = (-c0 * inv(c1)) % P
        return val

    def viol(val, cutoff=10**9):
        vi = 0
        for atom in matoms:
            s = 0
            for coef, m in atom:
                t = coef
                for x in m: t = (t * val[x]) % P
                s = (s + t) % P
            if s:
                vi += 1
                if vi >= cutoff: return vi
        return vi

    t1 = time.time()
    v0 = viol(evaluate(set()))
    print(f"baseline DAG-eval violated={v0}, eval time={time.time()-t1:.2f}s", flush=True)

    # WalkSAT / annealing
    best = v0; bestS = set()
    S = set(); curV = v0
    it = 0
    restarts = 0
    while time.time() - t0 < budget:
        it += 1
        b = random.choice(control)
        nS = set(S); nS.discard(b) if b in nS else nS.add(b)
        nv = viol(evaluate(nS), cutoff=curV + 30)
        if nv <= curV or random.random() < 0.05:
            S, curV = nS, nv
            if curV < best:
                best = curV; bestS = set(S)
                print(f"  it{it} NEW BEST viol={best} |S|={len(S)} ({time.time()-t0:.0f}s)", flush=True)
                if best == 0:
                    # verify over Z with full propagation
                    print("  candidate with 0 DAG-violations! verifying over Z...", flush=True)
                    json.dump(sorted(bestS), open(f'fastwalk_solution_bits_{SEED}.json', 'w'))
                    val = evaluate(S)
                    json.dump({f"x_{i}": val[i] for i in range(NVARS)}, open(f'cand_fastwalk_modp_{SEED}.json', 'w'))
                    return
        if it % 200 == 0:
            # random restart if stuck
            if random.random() < 0.3:
                S = set(random.sample(control, random.randint(0, 6))); curV = viol(evaluate(S)); restarts += 1
            print(f"  it{it} cur={curV} best={best} |S|={len(S)} restarts={restarts} ({time.time()-t0:.0f}s)", flush=True)
    print(f"WalkSAT done: best {best} bits {sorted(bestS)} after {it} iters", flush=True)
    json.dump({"best": best, "bits": sorted(bestS), "iters": it}, open(f'fastwalk_result_{SEED}.json', 'w'))

if __name__ == '__main__':
    main()
