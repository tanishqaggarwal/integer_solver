#!/usr/bin/env python3
"""Simulated annealing over subsets of the 255 control bits. State = set of bits
forced to 1; energy = number of violated atoms after full propagation. Accepts
non-improving moves with annealing probability to escape the rugged plateau."""
import json, time, math
from collections import deque
from propagate import load_atoms, atom_vars
from repair import ProvEngine, boolean_vars

NVARS = 38748
# deterministic PRNG (Date/random are fine here; this is a standalone job)
import random
random.seed(12345)

def main():
    import sys
    budget_s = int(sys.argv[1]) if len(sys.argv) > 1 else 5400
    atoms = load_atoms()
    bset = boolean_vars(atoms)
    base = ProvEngine(atoms); base.propagate()
    bval = list(base.val); bprov = list(base.prov); bdom = dict(base.domain)
    control = json.load(open('control_bits.json'))

    def energy(ones):
        eng = ProvEngine(atoms)
        eng.val = list(bval); eng.prov = list(bprov); eng.domain = dict(bdom)
        eng.wl = deque(); eng.inwl = [False] * len(atoms)
        for b in ones:
            if eng.val[b] is None: eng.assign(b, 1, ('f', ()))
        eng.propagate()
        for v in [b for b in bset if eng.val[b] is None]:
            if eng.val[v] is None: eng.assign(v, 0, ('f', ())); eng.propagate()
        for v in range(NVARS):
            if eng.val[v] is None: eng.assign(v, 0, ('f', ())); eng.propagate()
        val = [x if x is not None else 0 for x in eng.val]
        viol = 0
        for poly in atoms:
            s = 0
            for m, c in poly.items():
                t = c
                for var in m: t *= val[var]
                s += t
            if s != 0: viol += 1
        return viol, val

    t0 = time.time()
    S = set()
    curE, _ = energy(S)
    bestE = curE; bestS = set(S)
    print(f"start energy (violated atoms) = {curE}", flush=True)
    it = 0
    while time.time() - t0 < budget_s:
        it += 1
        T = max(0.05, 2.5 * (1 - (time.time() - t0) / budget_s))  # cooling
        b = random.choice(control)
        newS = set(S)
        if b in newS: newS.discard(b)
        else: newS.add(b)
        newE, val = energy(newS)
        dE = newE - curE
        if dE <= 0 or random.random() < math.exp(-dE / T):
            S = newS; curE = newE
            if curE < bestE:
                bestE = curE; bestS = set(S)
                print(f"  it{it} T{T:.2f} NEW BEST energy={bestE} |S|={len(S)} ({time.time()-t0:.0f}s)", flush=True)
                if bestE == 0:
                    json.dump({f"x_{i}": val[i] for i in range(NVARS)}, open('cand_sa_solved.json', 'w'))
                    json.dump(sorted(bestS), open('sa_solution_bits.json', 'w'))
                    print("SOLVED by SA! wrote cand_sa_solved.json", flush=True)
                    return
        if it % 25 == 0:
            print(f"  it{it} T{T:.2f} cur={curE} best={bestE} |S|={len(S)} ({time.time()-t0:.0f}s)", flush=True)
    print(f"SA finished: best energy {bestE} with {len(bestS)} bits after {it} iters", flush=True)
    json.dump({"bestE": bestE, "bestS": sorted(bestS), "iters": it}, open('sa_result.json', 'w'))

if __name__ == '__main__':
    main()
