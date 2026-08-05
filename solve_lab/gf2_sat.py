#!/usr/bin/env python3
"""Encode the mod-2 (GF(2)) reduction of the whole atom system as CNF and solve
with a SAT solver. Each var is a bit; monomials are AND-defined aux vars; each
atom is XOR(monomials)=const via a Tseitin chain. A core bit's mod-2 value is its
actual 0/1 value, so a SAT model gives candidate core bits (verify over Z)."""
import json, time, sys
from collections import defaultdict
from pysat.solvers import Cadical153 as SAT
from pysat.formula import IDPool

def monokey(m):
    return tuple(sorted(set(m)))

def main():
    t0 = time.time()
    atoms = [json.loads(l)['poly'] for l in open('atoms/poly_atoms.jsonl')]
    NV = 38748
    pool = IDPool()
    def xvar(i): return pool.id(('x', i))
    def mvar(mono): return pool.id(('m', mono))

    clauses = []
    # AND definitions for monomials
    mono_seen = set()
    def define_mono(mono):
        if mono in mono_seen: return
        mono_seen.add(mono)
        if len(mono) == 1:
            return  # monomial var IS the x var
        mv = mvar(mono)
        lits = [xvar(v) for v in mono]
        # mv <-> AND(lits)
        for l in lits:
            clauses.append([-mv, l])
        clauses.append([mv] + [-l for l in lits])

    def mono_lit(mono):
        if len(mono) == 1:
            return xvar(mono[0])
        return mvar(mono)

    unsat_const = 0
    n_atoms = 0
    for poly in atoms:
        acc = defaultdict(int)
        for m, c in poly:
            acc[monokey(tuple(m))] ^= (c & 1)
        terms = []; const = 0
        for m, b in acc.items():
            if not b: continue
            if m == ():
                const ^= 1
            else:
                define_mono(m)
                terms.append(mono_lit(m))
        if not terms:
            if const: unsat_const += 1  # 0 = 1 mod 2 -> infeasible
            continue
        n_atoms += 1
        # XOR(terms) = const  via Tseitin chain
        if len(terms) == 1:
            clauses.append([terms[0]] if const else [-terms[0]])
        else:
            prev = terms[0]
            for i in range(1, len(terms)):
                t = terms[i]
                if i < len(terms) - 1:
                    nxt = pool.id(('xor', n_atoms, i))
                else:
                    nxt = None
                if nxt is None:
                    # prev XOR t = const
                    if const:
                        clauses.append([prev, t]); clauses.append([-prev, -t])
                    else:
                        clauses.append([prev, -t]); clauses.append([-prev, t])
                else:
                    # nxt = prev XOR t
                    clauses.append([-nxt, -prev, -t]); clauses.append([-nxt, prev, t])
                    clauses.append([nxt, -prev, t]); clauses.append([nxt, prev, -t])
                    prev = nxt
    print(f"encoded: {n_atoms} atoms, {len(mono_seen)} monomials, {len(clauses)} clauses, "
          f"{pool.top} vars, unsat_const_rows={unsat_const} ({time.time()-t0:.1f}s)", flush=True)
    if unsat_const:
        print("NOTE: some atoms are 0=1 mod 2 -> the all-atoms-zero system is INFEASIBLE mod 2!")

    s = SAT(bootstrap_with=clauses)
    print("solving GF(2) SAT...", flush=True)
    r = s.solve()
    print(f"SAT result: {r} ({time.time()-t0:.0f}s)", flush=True)
    if r:
        model = set(l for l in s.get_model() if l > 0)
        mainbits = json.load(open('main_comp.json'))['main_bits']
        core = {b: (1 if xvar(b) in model else 0) for b in mainbits}
        from collections import Counter
        print("core bit parity distribution:", dict(Counter(core.values())))
        full = {f"x_{i}": (1 if xvar(i) in model else 0) for i in range(NV)}
        json.dump(full, open('gf2_sat_model.json', 'w'))
        json.dump({str(k): v for k, v in core.items()}, open('gf2_sat_corebits.json', 'w'))
        print("wrote gf2_sat_model.json, gf2_sat_corebits.json")

if __name__ == '__main__':
    main()
