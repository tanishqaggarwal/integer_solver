#!/usr/bin/env python3
"""Iterative unsat-core solve (implicit hitting set). Start with all free
value-inputs = 0 (UNSAT); the unsat core names value-inputs that must be nonzero.
Free those, re-solve; if still UNSAT, add the new core. Converges to the set of
nonzero value-inputs. Freed value-inputs are unbounded integers; bits are {0,1}."""
import json, time, sys, re
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
    free_vi = set(v for v in genuine if eng.val[v] is None)

    zvars = {}
    def zv(v):
        if v not in zvars: zvars[v] = z3.Int(f"x_{v}")
        return zvars[v]
    s = z3.Solver(); s.set("timeout", 600*1000)
    for b in mainbits:
        s.add(z3.Or(zv(b) == 0, zv(b) == 1))
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
        s.add(expr == 0)
    print(f"model built ({time.time()-t0:.0f}s), free_vi={len(free_vi)}", flush=True)

    F = set()   # value-inputs allowed nonzero
    for it in range(40):
        assumptions = [zv(v) == 0 for v in free_vi if v not in F]
        r = s.check(assumptions)
        print(f"iter {it}: |F|={len(F)} -> {r} ({time.time()-t0:.0f}s)", flush=True)
        if r == z3.sat:
            m = s.model(); out = dict(cand)
            for v in mainv:
                out[f"x_{v}"] = int(m.eval(zv(v), model_completion=True).as_long())
            json.dump(out, open('cand_z3iter.json', 'w'))
            nz = [v for v in free_vi if out[f"x_{v}"] != 0]
            print(f"*** SAT! nonzero value-inputs: {nz}. wrote cand_z3iter.json — VERIFY ***", flush=True)
            return
        if r != z3.unsat:
            print("unknown/timeout — stopping", flush=True); return
        core = s.unsat_core()
        newv = set()
        for a in core:
            mm = re.search(r'x_(\d+)', str(a))
            if mm: newv.add(int(mm.group(1)))
        add = newv - F
        if not add:
            print("no new core vars — stuck", flush=True); return
        F |= add
        print(f"  core adds {sorted(add)}; |F| now {len(F)}", flush=True)
        json.dump(sorted(F), open('z3iter_F.json', 'w'))

if __name__ == '__main__':
    main()
