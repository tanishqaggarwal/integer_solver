#!/usr/bin/env python3
"""z3 on the main component with free value-inputs constrained to a FINITE domain
{0} ∪ {a few residues} and bits to {0,1}. Finite domains can make the otherwise
unbounded-integer NIA problem decidable. Non-main vars fixed to the near-solution."""
import json, time, sys
import z3
from propagate import load_atoms
from repair import ProvEngine, boolean_vars

NVARS = 38748

def main():
    tout = int(sys.argv[1]) if len(sys.argv) > 1 else 2400
    t0 = time.time()
    atoms = [json.loads(l)['poly'] for l in open('atoms/poly_atoms.jsonl')]
    mc = json.load(open('main_comp.json')); mainv = set(mc['main_vars']); mainbits = set(mc['main_bits'])
    cand = json.load(open('best/best_partial_39013.json'))
    fixed = [0]*NVARS
    for k, v in cand.items(): fixed[int(k[2:])] = v

    # free value-inputs = genuine free (never clean target, not huge x_B, undetermined)
    summ = json.load(open('atoms/summary.json')); inputs = set(summ['inputs'])
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
    # determine which genuine are truly free (undetermined by pins+bits=0)
    eng = ProvEngine(load_atoms()); eng.propagate()
    bset0 = boolean_vars(load_atoms())
    for b in [x for x in bset0 if eng.val[x] is None]:
        if eng.val[b] is None: eng.assign(b, 0, ('f', ())); eng.propagate()
    free_vi = [v for v in genuine if eng.val[v] is None]

    # candidate residue domain (the distinct nonzero residue values seen in near-solution)
    residues = sorted(set(abs(v) for v in cand.values() if 200 <= abs(int(v)).bit_length() <= 300))
    print(f"free value-inputs: {len(free_vi)}, residue domain size {len(residues)+1} (0 + {len(residues)})", flush=True)

    zvars = {}
    def zv(v):
        if v not in zvars: zvars[v] = z3.Int(f"x_{v}")
        return zvars[v]

    s = z3.Solver(); s.set("timeout", tout*1000)
    # domain constraints
    for v in free_vi:
        s.add(z3.Or([zv(v) == 0] + [zv(v) == r for r in residues]))
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
    print(f"z3 model: {len(zvars)} vars, {ncon} constraints ({time.time()-t0:.0f}s); solving...", flush=True)
    r = s.check()
    print(f"z3 result: {r} ({time.time()-t0:.0f}s)", flush=True)
    if r == z3.sat:
        m = s.model(); out = dict(cand)
        for v in mainv:
            val = m.eval(zv(v), model_completion=True)
            out[f"x_{v}"] = int(val.as_long())
        json.dump(out, open('cand_z3residue.json', 'w'))
        print("*** SAT! wrote cand_z3residue.json — verify over Z ***", flush=True)

if __name__ == '__main__':
    main()
