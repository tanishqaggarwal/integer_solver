#!/usr/bin/env python3
"""Targeted finite-domain z3 on the fan-in cone of the 4 violated atoms.
Free: 178 value-inputs (domain {0}∪residues), 243 bits ({0,1}), 2859 wires (int).
Everything outside the cone is fixed to the near-solution. Require every atom that
touches a cone var to be 0. SAT => a full solution (verify over Z)."""
import json, time, sys
import z3

NVARS = 38748

def main():
    tout = int(sys.argv[1]) if len(sys.argv) > 1 else 1800
    t0 = time.time()
    atoms = [json.loads(l)['poly'] for l in open('atoms/poly_atoms.jsonl')]
    cone = json.load(open('cone4.json'))
    free_vi = set(cone['free_vi']); bits = set(cone['bits']); wires = set(cone['wires'])
    conevars = free_vi | bits | wires
    cand = json.load(open('best/best_partial_39013.json'))
    fixed = [0]*NVARS
    for k, v in cand.items(): fixed[int(k[2:])] = v
    # residue domain: all distinct HUGE constants + values seen in near-solution (~290-bit)
    residues = sorted(set(int(x) for x in open('big_literals.txt').read().split()))
    print(f"cone: {len(conevars)} vars ({len(free_vi)} value-inputs, {len(bits)} bits); residues {len(residues)}", flush=True)

    zvars = {}
    def zv(v):
        if v not in zvars: zvars[v] = z3.Int(f"x_{v}")
        return zvars[v]
    s = z3.Solver(); s.set("timeout", tout*1000)
    for v in free_vi:
        s.add(z3.Or([zv(v) == 0] + [zv(v) == r for r in residues] + [zv(v) == -r for r in residues]))
    for b in bits:
        s.add(z3.Or(zv(b) == 0, zv(b) == 1))
    ncon = 0
    for poly in atoms:
        vs = set()
        for m, c in poly: vs.update(m)
        if not (vs & conevars): continue
        expr = z3.IntVal(0)
        for m, c in poly:
            term = z3.IntVal(c)
            for var in m:
                term = term * (zv(var) if var in conevars else z3.IntVal(fixed[var]))
            expr = expr + term
        s.add(expr == 0); ncon += 1
    print(f"z3: {len(zvars)} vars, {ncon} constraints built ({time.time()-t0:.0f}s); solving...", flush=True)
    r = s.check()
    print(f"z3 result: {r} ({time.time()-t0:.0f}s)", flush=True)
    if r == z3.sat:
        m = s.model(); out = dict(cand)
        for v in conevars:
            val = m.eval(zv(v), model_completion=True)
            out[f"x_{v}"] = int(val.as_long())
        json.dump(out, open('cand_z3cone.json', 'w'))
        print("*** SAT! wrote cand_z3cone.json — VERIFY over Z with checker.py ***", flush=True)
    elif r == z3.unsat:
        print("UNSAT: witness changes vars outside this cone, or value-inputs aren't residue-or-0", flush=True)

if __name__ == '__main__':
    main()
