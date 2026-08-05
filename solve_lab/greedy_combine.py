#!/usr/bin/env python3
"""Greedy multi-bit combiner: from all-zero-bits baseline, greedily add bit
flips (=1) that most reduce the number of violated atoms, testing actual
combinations (bit effects interact through the circuit)."""
import json, time
from collections import deque
from propagate import load_atoms, atom_vars, substitute, solve_single
from repair import ProvEngine, boolean_vars

NVARS = 38748

def main():
    t0 = time.time()
    atoms = load_atoms()
    bset = boolean_vars(atoms)
    base = ProvEngine(atoms); base.propagate()
    base_val = list(base.val); base_prov = list(base.prov); base_dom = dict(base.domain)

    def run(first_ones):
        eng = ProvEngine(atoms)
        eng.val = list(base_val); eng.prov = list(base_prov); eng.domain = dict(base_dom)
        eng.wl = deque(); eng.inwl = [False] * len(atoms)
        for b in first_ones:
            if eng.val[b] is None: eng.assign(b, 1, ('flip', ()))
        eng.propagate()
        for v in [b for b in bset if eng.val[b] is None]:
            if eng.val[v] is None: eng.assign(v, 0, ('free', ())); eng.propagate()
        for v in range(NVARS):
            if eng.val[v] is None: eng.assign(v, 0, ('free', ())); eng.propagate()
        val = [x if x is not None else 0 for x in eng.val]
        viol = 0; vset = []
        for ai, poly in enumerate(atoms):
            s = 0
            for m, c in poly.items():
                t = c
                for var in m: t *= val[var]
                s += t
            if s != 0:
                viol += 1; vset.append(ai)
        return viol, vset, val

    # candidate pool: improving single flips
    res = json.load(open('flip_results.json'))
    pool = [b for (viol, b, nc) in res if viol < 4]
    print(f"candidate pool: {len(pool)} improving bits; baseline viol=4")

    chosen = []
    cur_viol, cur_vset, cur_val = run(chosen)
    print(f"start viol={cur_viol} atoms={cur_vset}")
    for step in range(8):
        best = None
        for b in pool:
            if b in chosen: continue
            viol, vset, val = run(chosen + [b])
            if best is None or viol < best[0]:
                best = (viol, b, vset, val)
        if best is None or best[0] >= cur_viol:
            print(f"step {step}: no improvement (best {best[0] if best else '-'}); stop")
            break
        cur_viol, b, cur_vset, cur_val = best
        chosen.append(b)
        print(f"step {step}: +x_{b}=1 -> viol={cur_viol} atoms={cur_vset}  ({time.time()-t0:.0f}s)", flush=True)
        if cur_viol == 0:
            print("SOLVED via greedy flips!")
            json.dump({f"x_{i}": cur_val[i] for i in range(NVARS)}, open('cand_greedy_solved.json', 'w'))
            break
    print(f"chosen bits: {chosen}, final viol {cur_viol}")
    json.dump({"chosen": chosen, "viol": cur_viol, "vset": cur_vset}, open('greedy_result.json', 'w'))

if __name__ == '__main__':
    main()
