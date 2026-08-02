#!/usr/bin/env python3
"""Build the GF(2) CNF once, then for each core bit test whether it is FORCED
mod 2 (i.e. the opposite value is UNSAT). Forced-mod-2 bits are definite witness
bit values (the Z-witness satisfies the mod-2 system)."""
import json, time, sys
from collections import defaultdict
from pysat.solvers import Cadical153 as SAT
from pysat.formula import IDPool

def monokey(m): return tuple(sorted(set(m)))

def build():
    atoms = [json.loads(l)['poly'] for l in open('atoms/poly_atoms.jsonl')]
    pool = IDPool()
    def xvar(i): return pool.id(('x', i))
    def mvar(mono): return pool.id(('m', mono))
    clauses = []; mono_seen = set()
    def define_mono(mono):
        if mono in mono_seen or len(mono) == 1:
            mono_seen.add(mono); return
        mono_seen.add(mono); mv = mvar(mono); lits = [xvar(v) for v in mono]
        for l in lits: clauses.append([-mv, l])
        clauses.append([mv] + [-l for l in lits])
    def mono_lit(mono): return xvar(mono[0]) if len(mono) == 1 else mvar(mono)
    n = 0
    for poly in atoms:
        acc = defaultdict(int)
        for m, c in poly: acc[monokey(tuple(m))] ^= (c & 1)
        terms = []; const = 0
        for m, b in acc.items():
            if not b: continue
            if m == (): const ^= 1
            else: define_mono(m); terms.append(mono_lit(m))
        if not terms: continue
        n += 1
        if len(terms) == 1:
            clauses.append([terms[0]] if const else [-terms[0]])
        else:
            prev = terms[0]
            for i in range(1, len(terms)):
                t = terms[i]
                nxt = pool.id(('xor', n, i)) if i < len(terms) - 1 else None
                if nxt is None:
                    if const: clauses += [[prev, t], [-prev, -t]]
                    else: clauses += [[prev, -t], [-prev, t]]
                else:
                    clauses += [[-nxt,-prev,-t],[-nxt,prev,t],[nxt,-prev,t],[nxt,prev,-t]]
                    prev = nxt
    return clauses, xvar

def main():
    t0 = time.time()
    clauses, xvar = build()
    mainbits = json.load(open('main_comp.json'))['main_bits']
    s = SAT(bootstrap_with=clauses)
    assert s.solve(), "base GF(2) UNSAT?!"
    model = set(l for l in s.get_model() if l > 0)
    print(f"CNF built + base solve in {time.time()-t0:.0f}s", flush=True)
    forced = {}
    for k, b in enumerate(mainbits):
        cur = 1 if xvar(b) in model else 0
        # try opposite
        assume = -xvar(b) if cur == 1 else xvar(b)
        if not s.solve(assumptions=[assume]):
            forced[b] = cur   # opposite is UNSAT -> forced to cur
        else:
            # refresh model occasionally
            model = set(l for l in s.get_model() if l > 0)
        if (k+1) % 32 == 0:
            print(f"   tested {k+1}/256, forced so far {len(forced)} ({time.time()-t0:.0f}s)", flush=True)
    from collections import Counter
    print(f"\nCORE BITS FORCED MOD 2: {len(forced)}/256, values {dict(Counter(forced.values()))}")
    json.dump({str(k): v for k, v in forced.items()}, open('gf2_forced_core.json', 'w'))
    print(f"free-mod-2 core bits: {256-len(forced)}  ({time.time()-t0:.0f}s)")

if __name__ == '__main__':
    main()
