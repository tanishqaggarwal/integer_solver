#!/usr/bin/env python3
"""Solve the main residual component with z3. Non-main variables are fixed to
their near-solution (cand_forward2) values; main-component variables are z3
Ints. Every atom touching a main var becomes a == 0 constraint."""
import json, time, sys
import z3

NVARS = 38748
TIMEOUT_MS = int(sys.argv[1]) * 1000 if len(sys.argv) > 1 else 900_000

def main():
    t0 = time.time()
    atoms = [json.loads(l)['poly'] for l in open('solve_lab/atoms/poly_atoms.jsonl')]
    mc = json.load(open('solve_lab/main_comp.json'))
    mainv = set(mc['main_vars'])
    mainbits = set(mc['main_bits'])
    cand = json.load(open('solve_lab/cand_forward2.json'))
    fixed = [0] * NVARS
    for k, x in cand.items():
        fixed[int(k[2:])] = x

    zvars = {}
    def zv(v):
        if v not in zvars:
            zvars[v] = z3.Int(f"x_{v}")
        return zvars[v]

    s = z3.Solver()
    s.set("timeout", TIMEOUT_MS)

    n_constraints = 0
    for poly in atoms:
        vs = set()
        for m, c in poly: vs.update(m)
        if not (vs & mainv):
            continue  # atom fully outside main component -> already satisfied
        # build z3 expression, substituting fixed values for non-main vars
        expr = z3.IntVal(0)
        for m, c in poly:
            term = z3.IntVal(c)
            for var in m:
                term = term * (zv(var) if var in mainv else z3.IntVal(fixed[var]))
            expr = expr + term
        s.add(expr == 0)
        n_constraints += 1

    for b in mainbits:
        s.add(zv(b) >= 0, zv(b) <= 1)

    print(f"built z3 model: {len(zvars)} vars, {n_constraints} constraints in {time.time()-t0:.1f}s", flush=True)
    print(f"solving with timeout {TIMEOUT_MS/1000:.0f}s ...", flush=True)
    r = s.check()
    print(f"z3 result: {r}  ({time.time()-t0:.1f}s)", flush=True)
    if r == z3.sat:
        model = s.model()
        out = dict(cand)
        for v in mainv:
            val = model.eval(zv(v), model_completion=True)
            out[f"x_{v}"] = int(val.as_long()) if val is not None else fixed[v]
        json.dump(out, open('solve_lab/cand_z3main.json', 'w'))
        print("wrote cand_z3main.json", flush=True)
    else:
        print("no model (unsat/unknown/timeout)", flush=True)

if __name__ == '__main__':
    main()
