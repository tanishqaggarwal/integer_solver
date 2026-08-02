#!/usr/bin/env python3
"""Fast topological forward-evaluator for the whole system.

Each wire is computed once, in dependency order, from its defining atom.
Huge atoms bit*(x_B-HUGE)=s*x_C define x_B conditionally: x_B = HUGE+s*x_C if
bit=1 else x_B is a free input (0). Given a bit assignment this evaluates all
wires with one pass, then counts violated atoms. Much faster than re-propagation.
"""
import json, time, sys
from collections import defaultdict, deque
from propagate import load_atoms, atom_vars, substitute, solve_single
from repair import ProvEngine, boolean_vars

NVARS = 38748

def build():
    atoms = load_atoms()
    bset = boolean_vars(atoms)
    # provenance at bits=0
    eng = ProvEngine(atoms); eng.propagate()
    base_forced = {i: eng.val[i] for i in range(NVARS) if eng.val[i] is not None}
    # continue: zero free bits, zero-fill, to get full provenance DAG for non-huge wires
    prov = list(eng.prov)  # snapshot after pins
    # We rebuild a clean DAG below using ALL definitional atoms.

    # identify huge atoms and their (bit, xB, defining poly)
    huge_def = {}   # xB -> (atom_index, bit, poly)
    huge_atoms = set()
    for ai, poly in enumerate(atoms):
        d = dict(poly)
        # find degree-1 monomial with huge coef -> that var is bit, coef is -HUGE
        bigs = [(m, c) for m, c in poly.items() if len(m) == 1 and abs(c) >= 10**20]
        deg2 = [(m, c) for m, c in poly.items() if len(m) == 2]
        if len(bigs) == 1 and deg2:
            bit = bigs[0][0][0]
            # x_B: the var paired with bit in a deg2 monomial coef +-1
            xB = None
            for m, c in deg2:
                if bit in m and abs(c) == 1:
                    xB = m[0] if m[1] == bit else m[1]
                    break
            if xB is not None:
                huge_def[xB] = (ai, bit, poly)
                huge_atoms.add(ai)
    print(f"huge-atom loads: {len(huge_def)} x_B wires, over {len(huge_atoms)} atoms", file=sys.stderr)

    # definition of each wire: from provenance of a FULL solve (bits=0), except x_B via huge
    eng2 = ProvEngine(atoms); eng2.propagate()
    bits_all = [b for b in bset if eng2.val[b] is None]
    for b in bits_all:
        if eng2.val[b] is None: eng2.assign(b, 0, ('free', ())); eng2.propagate()
    for v in range(NVARS):
        if eng2.val[v] is None: eng2.assign(v, 0, ('free', ())); eng2.propagate()

    define = {}   # wire -> ('atom', ai, deps) or ('huge', ai, bit) or ('input',)
    for v in range(NVARS):
        p = eng2.prov[v]
        if v in huge_def:
            ai, bit, poly = huge_def[v]
            deps = tuple(x for x in atom_vars(atoms[ai]) if x != v)
            define[v] = ('huge', ai, deps)
        elif p is None:
            define[v] = ('input',)
        elif p[0] == 'free':
            define[v] = ('input',)
        else:
            ai, deps = p
            define[v] = ('atom', ai, deps)

    # topological order
    color = [0] * NVARS  # 0=unvisited,1=in progress,2=done
    order = []
    def visit_iter(start):
        stack = [(start, False)]
        while stack:
            v, processed = stack.pop()
            if processed:
                if color[v] != 2:
                    color[v] = 2; order.append(v)
                continue
            if color[v] == 2: continue
            if color[v] == 1:  # cycle - break by treating as input
                continue
            color[v] = 1
            stack.append((v, True))
            kind = define[v]
            if kind[0] in ('atom', 'huge'):
                for d in kind[2] if kind[0] == 'huge' else kind[2]:
                    if color[d] == 0:
                        stack.append((d, False))
    for v in range(NVARS):
        if color[v] == 0:
            visit_iter(v)

    return atoms, bset, define, order, base_forced, eng2.val

def make_evaluator(atoms, define, order):
    codes = {}
    for v in range(NVARS):
        k = define[v]
        if k[0] == 'atom':
            codes[v] = ('atom', atoms[k[1]], v)
        elif k[0] == 'huge':
            codes[v] = ('huge', atoms[k[1]], v)
    return codes

def main():
    t0 = time.time()
    atoms, bset, define, order, base_forced, ref_val = build()
    print(f"built DAG + topo order ({len(order)} wires) in {time.time()-t0:.1f}s", file=sys.stderr)

    mainbits = json.load(open('main_comp.json'))['main_bits']

    def evaluate(one_bits):
        val = [None] * NVARS
        oneset = set(one_bits)
        for v in order:
            k = define[v]
            if k[0] == 'input':
                if v in bset:
                    val[v] = 1 if v in oneset else 0
                else:
                    val[v] = 0
                continue
            poly = atoms[k[1]]
            # substitute known deps, solve for v
            red = {}
            for m, c in poly.items():
                coef = c; newm = []
                ok = True
                for x in m:
                    if x == v:
                        newm.append(x)
                    elif val[x] is None:
                        ok = False; break
                    else:
                        coef *= val[x]
                if not ok:
                    red = None; break
                key = tuple(newm)
                red[key] = red.get(key, 0) + coef
            if red is None:
                val[v] = 0; continue
            # red is c1*v + c0 (+ maybe c2*v^2). solve linear
            c0 = red.get((), 0); c1 = red.get((v,), 0); c2 = red.get((v, v), 0)
            if c2 == 0 and c1 != 0 and (-c0) % c1 == 0:
                val[v] = (-c0) // c1
            elif c2 == 0 and c1 != 0:
                val[v] = 0  # non-integer -> treat as 0 (inconsistent)
            else:
                val[v] = 0
        return val

    def count_viol(val):
        viol = 0
        for poly in atoms:
            s = 0
            for m, c in poly.items():
                t = c
                for x in m: t *= val[x]
                s += t
            if s != 0: viol += 1
        return viol

    # verify against reference (bits=0 -> x_24550 forced to 1 in ref; here forced by eval)
    t1 = time.time()
    v0 = evaluate([])
    viol0 = count_viol(v0)
    # compare to ref on forced vars
    match = sum(1 for i in range(NVARS) if v0[i] == (ref_val[i] if ref_val[i] is not None else 0))
    print(f"eval bits=0: {viol0} violated atoms, {match}/{NVARS} match near-solution, eval={time.time()-t1:.2f}s")

    if len(sys.argv) > 1 and sys.argv[1] == 'search':
        import itertools
        control = json.load(open('control_bits.json'))
        best = (viol0, [])
        # pairs over improving bits
        res = json.load(open('flip_results.json'))
        pool = [b for (vv, b, nc) in res if vv <= 4][:60]
        print(f"searching pairs over {len(pool)} bits...")
        n = 0
        for a, b in itertools.combinations(pool, 2):
            val = evaluate([a, b]); vc = count_viol(val); n += 1
            if vc < best[0]:
                best = (vc, [a, b]); print(f"  pair x_{a},x_{b}: {vc} violated", flush=True)
                if vc == 0:
                    json.dump({f"x_{i}": val[i] for i in range(NVARS)}, open('cand_fast_solved.json', 'w'))
                    print("SOLVED"); return
        print(f"tested {n} pairs, best {best}")

if __name__ == '__main__':
    main()
