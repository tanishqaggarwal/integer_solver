#!/usr/bin/env python3
"""Accurate fast evaluator: replay propagation's recorded assignment order in a
single mod-P pass (valid feed-forward, matches propagation -> baseline ~4-5, and
0 exactly at a true solution). Then WalkSAT/anneal over the core bits."""
import json, time, sys, random
from fast_eval import NVARS  # just for NVARS
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

def main():
    budget = int(sys.argv[1]) if len(sys.argv) > 1 else 6000
    SEED = int(sys.argv[2]) if len(sys.argv) > 2 else 1
    random.seed(SEED)
    t0 = time.time()
    atoms = load_atoms(); bset = boolean_vars(atoms)
    eo = json.load(open('eval_order.json'))
    order = eo['order']; prov = eo['prov']
    mainv = set(json.load(open('main_comp.json'))['main_vars'])
    control = json.load(open('control_bits.json'))

    defatom = {}
    for v in range(NVARS):
        p = prov[v]
        if p and p[0] >= 0:
            defatom[v] = [(c % P, tuple(m)) for m, c in atoms[p[0]].items()]
    matoms = [[(c % P, tuple(m)) for m, c in poly.items()]
              for poly in atoms if set().union(*[set(m) for m in poly]) & mainv]
    is_bit = [False] * NVARS
    for b in bset: is_bit[b] = True
    print(f"loaded; order {len(order)}, main atoms {len(matoms)} ({time.time()-t0:.1f}s)", flush=True)

    def evaluate(oneset):
        val = [0] * NVARS
        for v in order:
            da = defatom.get(v)
            if da is None:
                val[v] = 1 if (is_bit[v] and v in oneset) else 0
                continue
            c0 = c1 = 0
            for coef, m in da:
                if v in m:
                    t = coef
                    for x in m:
                        if x != v: t = (t * val[x]) % P
                    c1 = (c1 + t) % P
                else:
                    t = coef
                    for x in m: t = (t * val[x]) % P
                    c0 = (c0 + t) % P
            if c1: val[v] = (-c0 * inv(c1)) % P
            else: val[v] = 0
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
    print(f"baseline (bits=0) accurate DAG viol = {v0}  (Z propagation gives 4), eval {time.time()-t1:.2f}s", flush=True)

    best = v0; bestS = set(); S = set(); curV = v0; it = 0; restarts = 0
    while time.time() - t0 < budget:
        it += 1
        b = random.choice(control)
        nS = set(S)
        nS.discard(b) if b in nS else nS.add(b)
        nv = viol(evaluate(nS), cutoff=curV + 20)
        if nv <= curV or random.random() < 0.04:
            S, curV = nS, nv
            if curV < best:
                best = curV; bestS = set(S)
                print(f"  it{it} NEW BEST viol={best} |S|={len(S)} ({time.time()-t0:.0f}s)", flush=True)
                if best == 0:
                    json.dump(sorted(bestS), open(f'fw2_solution_bits_{SEED}.json', 'w'))
                    val = evaluate(S)
                    json.dump({f"x_{i}": val[i] for i in range(NVARS)}, open(f'cand_fw2_{SEED}.json', 'w'))
                    print("  *** 0 DAG-violations candidate saved -> verify over Z ***", flush=True)
                    return
        if it % 300 == 0:
            if random.random() < 0.35:
                S = set(random.sample(control, random.randint(0, 8)))
                curV = viol(evaluate(S)); restarts += 1
            print(f"  it{it} cur={curV} best={best} |S|={len(S)} r={restarts} ({time.time()-t0:.0f}s)", flush=True)
    print(f"done best={best} bits={sorted(bestS)} iters={it}", flush=True)
    json.dump({"best": best, "bits": sorted(bestS), "iters": it}, open(f'fw2_result_{SEED}.json', 'w'))

if __name__ == '__main__':
    main()
