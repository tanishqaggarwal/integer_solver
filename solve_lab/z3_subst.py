#!/usr/bin/env python3
"""Substituted-model solve: free value-inputs = F (residue domain); ALL other free
value-inputs are substituted as 0 (collapsing their products -> much smaller, faster
model). Bits free. If SAT verify; if UNSAT, grow F via an assumption-model unsat core.
Iterates."""
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
    residues = sorted(set(int(x) for x in open('big_literals.txt').read().split()))
    print(f"free_vi={len(free_vi)}, residues={len(residues)}", flush=True)

    F = set(json.load(open('unsat_core_vinputs.json')))  # start from first core
    print(f"initial F (from core) = {sorted(F)}", flush=True)

    for it in range(60):
        # substituted model: non-F free value-inputs -> 0
        zero_vi = free_vi - F
        s = z3.Solver(); s.set("timeout", 900*1000)
        zvars = {}
        def zv(v):
            if v not in zvars: zvars[v] = z3.Int(f"x_{v}")
            return zvars[v]
        def term_var(var):
            if var in zero_vi: return z3.IntVal(0)
            if var in mainv: return zv(var)
            return z3.IntVal(fixed[var])
        for v in F:
            s.add(z3.Or([zv(v) == 0] + [zv(v) == r for r in residues]))
        for b in mainbits:
            s.add(z3.Or(zv(b) == 0, zv(b) == 1))
        # tracked "F-member == 0" so we can grow via core if UNSAT
        track = {}
        for v in F:
            p = z3.Bool(f"t_{v}"); track[p] = v
            s.add(z3.Implies(p, zv(v) == 0))
        ncon = 0
        for poly in atoms:
            vs = set()
            for m, c in poly: vs.update(m)
            if not (vs & mainv): continue
            if vs & zero_vi and not (vs & (F | mainbits)) and not (vs - zero_vi - {x for x in vs if x not in mainv}):
                pass
            expr = z3.IntVal(0)
            for m, c in poly:
                t = z3.IntVal(c)
                for var in m: t = t * term_var(var)
                expr = expr + t
            s.add(expr == 0); ncon += 1
        r = s.check()
        print(f"iter {it}: |F|={len(F)} zero={len(zero_vi)} cons={ncon} -> {r} ({time.time()-t0:.0f}s)", flush=True)
        if r == z3.sat:
            m = s.model(); out = dict(cand)
            for v in mainv:
                vv = zvars.get(v)
                out[f"x_{v}"] = int(m.eval(vv, model_completion=True).as_long()) if vv is not None else (0 if v in zero_vi else fixed[v])
            json.dump(out, open('cand_z3subst.json', 'w'))
            nz = [v for v in F if out[f"x_{v}"] != 0]
            print(f"*** SAT! nonzero core value-inputs: {nz} — wrote cand_z3subst.json — VERIFY ***", flush=True)
            return
        if r != z3.unsat:
            print("unknown/timeout — stopping", flush=True); return
        # UNSAT: which zero_vi must be nonzero? re-derive core over an assumption model
        # cheaper: add all zero_vi as candidates via a fresh core call
        s2 = z3.Solver(); s2.set("timeout", 300*1000)
        zv2 = {}
        def z2(v):
            if v not in zv2: zv2[v] = z3.Int(f"x_{v}")
            return zv2[v]
        for b in mainbits: s2.add(z3.Or(z2(b) == 0, z2(b) == 1))
        for v in F:
            s2.add(z3.Or([z2(v) == 0] + [z2(v) == r for r in residues]))
        for poly in atoms:
            vs = set()
            for m, c in poly: vs.update(m)
            if not (vs & mainv): continue
            expr = z3.IntVal(0)
            for m, c in poly:
                t = z3.IntVal(c)
                for var in m: t = t * (z2(var) if var in mainv else z3.IntVal(fixed[var]))
                expr = expr + t
            s2.add(expr == 0)
        assum = [z2(v) == 0 for v in (free_vi - F)]
        r2 = s2.check(assum)
        if r2 != z3.unsat:
            print(f"grow-core check returned {r2} — stopping", flush=True); return
        core = s2.unsat_core()
        add = set()
        for a in core:
            mm = re.search(r'x_(\d+)', str(a))
            if mm: add.add(int(mm.group(1)))
        add -= F
        if not add:
            print("no new core — stopping", flush=True); return
        F |= add
        print(f"  grew F by {sorted(add)} -> |F|={len(F)}", flush=True)
        json.dump(sorted(F), open('z3subst_F.json', 'w'))

if __name__ == '__main__':
    main()
