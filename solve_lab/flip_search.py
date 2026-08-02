#!/usr/bin/env python3
"""Single-bit-flip search: from the all-zero-bits baseline, try setting each
control bit to 1 (first), re-propagate, and count violated atoms. Report bits
that reduce the violation count."""
import json, time, sys
from copy import copy
from collections import deque
from propagate import load_atoms, atom_vars, substitute, solve_single
from repair import ProvEngine, boolean_vars

NVARS = 38748

def build_evaluator(atoms):
    def ev_all(val):
        v = [x if x is not None else 0 for x in val]
        cnt = 0
        for poly in atoms:
            s = 0
            for m, c in poly.items():
                t = c
                for var in m: t *= v[var]
                s += t
            if s != 0: cnt += 1
        return cnt
    return ev_all

def run(atoms, base_val, base_prov, base_dom, bset, first_ones):
    eng = ProvEngine(atoms)
    eng.val = list(base_val)
    eng.prov = list(base_prov)
    eng.domain = dict(base_dom)
    eng.wl = deque()
    eng.inwl = [False] * len(atoms)
    # apply first_ones
    for b in first_ones:
        if eng.val[b] is None:
            eng.assign(b, 1, ('flip', ()))
    eng.propagate()
    for v in [b for b in bset if eng.val[b] is None]:
        if eng.val[v] is None:
            eng.assign(v, 0, ('free', ())); eng.propagate()
    for v in range(NVARS):
        if eng.val[v] is None:
            eng.assign(v, 0, ('free', ())); eng.propagate()
    return eng

def main():
    t0 = time.time()
    atoms = load_atoms()
    bset = boolean_vars(atoms)
    ev_all = build_evaluator(atoms)
    # baseline: initial propagation only
    base = ProvEngine(atoms)
    base.propagate()
    base_val = list(base.val); base_prov = list(base.prov); base_dom = dict(base.domain)
    print(f"baseline propagation done ({time.time()-t0:.1f}s), forced {sum(1 for x in base_val if x is not None)}")

    # baseline violations (all free bits -> 0)
    eng0 = run(atoms, base_val, base_prov, base_dom, bset, [])
    base_viol = ev_all(eng0.val)
    print(f"baseline violated atoms (all bits 0): {base_viol}")

    control = json.load(open('control_bits.json'))
    print(f"testing {len(control)} control-bit single flips...")
    results = []
    for i, b in enumerate(control):
        eng = run(atoms, base_val, base_prov, base_dom, bset, [b])
        viol = ev_all(eng.val)
        results.append((viol, b, len(eng.contra)))
        if viol < base_viol:
            print(f"  [IMPROVE] flip x_{b}=1 -> violated {viol} (contra {len(eng.contra)})", flush=True)
        if (i+1) % 40 == 0:
            print(f"   ...{i+1}/{len(control)} tested, {time.time()-t0:.0f}s", flush=True)
    results.sort()
    print("\nbest 15 single flips (violated_atoms, bit, contra):")
    for viol, b, nc in results[:15]:
        print(f"   x_{b}=1 -> {viol} violated, {nc} contra")
    json.dump(results, open('flip_results.json', 'w'))
    print(f"done {time.time()-t0:.0f}s; wrote flip_results.json")

if __name__ == '__main__':
    main()
