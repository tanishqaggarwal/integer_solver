#!/usr/bin/env python3
"""cvc5 (different SMT engine) on the localized substituted model: free the 4
core value-inputs (residue domain) + bits, substitute other free value-inputs=0.
cvc5's nonlinear-integer engine may synthesize a model where z3 stalls."""
import json, time, sys
import cvc5.pythonic as z3
from propagate import load_atoms
from repair import ProvEngine, boolean_vars
NVARS = 38748

def main():
    t0 = time.time()
    atoms = [json.loads(l)['poly'] for l in open('atoms/poly_atoms.jsonl')]
    mc = json.load(open('main_comp.json')); mainv = set(mc['main_vars']); mainbits = set(mc['main_bits'])
    cand = json.load(open('best/best_partial_39013.json')); fixed = [0]*NVARS
    for k, v in cand.items(): fixed[int(k[2:])] = v
    summ = json.load(open('atoms/summary.json')); inputs = set(summ['inputs'])
    A = load_atoms(); huge_xB = set()
    for poly in atoms:
        d = {tuple(m): c for m, c in poly}; bigs = [(m, c) for m, c in d.items() if len(m) == 1 and abs(c) >= 10**20]; deg2 = [(m, c) for m, c in d.items() if len(m) == 2]
        if len(bigs) == 1 and deg2:
            bit = bigs[0][0][0]
            for m, c in deg2:
                if bit in m and abs(c) == 1: huge_xB.add(m[0] if m[1] == bit else m[1]); break
    genuine = (inputs & mainv) - mainbits - huge_xB
    eng = ProvEngine(A); eng.propagate(); bset = boolean_vars(A)
    for b in [x for x in bset if eng.val[x] is None]:
        if eng.val[b] is None: eng.assign(b, 0, ('f', ())); eng.propagate()
    free_vi = set(v for v in genuine if eng.val[v] is None)
    residues = sorted(set(int(x) for x in open('big_literals.txt').read().split()))
    F = set(json.load(open('unsat_core_vinputs.json')))
    zero_vi = free_vi - F
    zvars = {}
    def zv(v):
        if v not in zvars: zvars[v] = z3.Int(f"x_{v}")
        return zvars[v]
    s = z3.Solver()
    for v in F: s.add(z3.Or([zv(v) == 0] + [zv(v) == r for r in residues]))
    for b in mainbits: s.add(z3.Or(zv(b) == 0, zv(b) == 1))
    for poly in atoms:
        vs = set()
        for m, c in poly: vs.update(m)
        if not (vs & mainv): continue
        expr = z3.IntVal(0)
        for m, c in poly:
            t = z3.IntVal(c)
            for var in m:
                t = t * (z3.IntVal(0) if var in zero_vi else (zv(var) if var in mainv else z3.IntVal(fixed[var])))
            expr = expr + t
        s.add(expr == 0)
    print(f"cvc5 model built ({time.time()-t0:.0f}s); solving...", flush=True)
    r = s.check()
    print(f"cvc5 result: {r} ({time.time()-t0:.0f}s)", flush=True)
    if str(r) == 'sat':
        m = s.model(); out = dict(cand)
        for v in mainv:
            try: out[f"x_{v}"] = int(str(m[zv(v)]))
            except Exception: pass
        json.dump(out, open('cand_cvc5.json', 'w'))
        print("*** cvc5 SAT! wrote cand_cvc5.json — VERIFY ***", flush=True)

if __name__ == '__main__':
    main()
