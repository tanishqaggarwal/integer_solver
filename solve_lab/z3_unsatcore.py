#!/usr/bin/env python3
"""Extract the UNSAT core of 'all free value-inputs = 0' (with bits free).
The core = the minimal subset of value-inputs that CANNOT all be 0 => exactly the
ones the witness must set nonzero. If small, we free only those (residue domain)."""
import json, time, sys
import z3
from propagate import load_atoms
from repair import ProvEngine, boolean_vars

NVARS = 38748

def main():
    t0 = time.time()
    atoms = [json.loads(l)['poly'] for l in open('atoms/poly_atoms.jsonl')]
    mc = json.load(open('main_comp.json')); mainv = set(mc['main_vars']); mainbits = set(mc['main_bits'])
    cand = json.load(open('best/best_partial_39013.json'))
    fixed = [0]*NVARS
    for k, v in cand.items(): fixed[int(k[2:])] = v
    summ = json.load(open('atoms/summary.json')); inputs = set(summ['inputs'])
    A = load_atoms()
    huge_xB = set()
    for poly in atoms:
        d = {tuple(m): c for m, c in poly}
        bigs = [(m, c) for m, c in d.items() if len(m) == 1 and abs(c) >= 10**20]
        deg2 = [(m, c) for m, c in d.items() if len(m) == 2]
        if len(bigs) == 1 and deg2:
            bit = bigs[0][0][0]
            for m, c in deg2:
                if bit in m and abs(c) == 1: huge_xB.add(m[0] if m[1] == bit else m[1]); break
    genuine = (inputs & mainv) - mainbits - huge_xB
    eng = ProvEngine(A); eng.propagate()
    bset = boolean_vars(A)
    for b in [x for x in bset if eng.val[x] is None]:
        if eng.val[b] is None: eng.assign(b, 0, ('f', ())); eng.propagate()
    free_vi = [v for v in genuine if eng.val[v] is None]
    print(f"free value-inputs: {len(free_vi)}", flush=True)

    zvars = {}
    def zv(v):
        if v not in zvars: zvars[v] = z3.Int(f"x_{v}")
        return zvars[v]
    s = z3.Solver(); s.set("timeout", 1200*1000)
    for b in mainbits:
        s.add(z3.Or(zv(b) == 0, zv(b) == 1))
    ncon = 0
    for poly in atoms:
        vs = set()
        for m, c in poly: vs.update(m)
        if not (vs & mainv): continue
        expr = z3.IntVal(0)
        for m, c in poly:
            term = z3.IntVal(c)
            for var in m:
                term = term * (zv(var) if var in mainv else z3.IntVal(fixed[var]))
            expr = expr + term
        s.add(expr == 0); ncon += 1
    # tracked assumptions: value-input == 0
    assumptions = [zv(v) == 0 for v in free_vi]
    print(f"model {len(zvars)} vars {ncon} cons built ({time.time()-t0:.0f}s); solving with assumptions...", flush=True)
    r = s.check(assumptions)
    print(f"result: {r} ({time.time()-t0:.0f}s)", flush=True)
    if r == z3.unsat:
        core = s.unsat_core()
        core_vars = []
        for a in core:
            sa = str(a)
            # a is "x_V == 0"
            import re
            m = re.search(r'x_(\d+)', sa)
            if m: core_vars.append(int(m.group(1)))
        print(f"UNSAT CORE: {len(core)} assumptions -> value-inputs that must be nonzero: {sorted(set(core_vars))}", flush=True)
        json.dump(sorted(set(core_vars)), open('unsat_core_vinputs.json', 'w'))
        print("wrote unsat_core_vinputs.json")
    elif r == z3.sat:
        print("SAT?! value-inputs=0 has a solution — unexpected", flush=True)

if __name__ == '__main__':
    main()
