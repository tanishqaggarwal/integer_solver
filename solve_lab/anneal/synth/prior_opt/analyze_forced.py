#!/usr/bin/env python3
"""
analyze_forced.py -- INDEPENDENT, from-scratch reproduction of the pre-QUBO
structural shrink, writing ALL outputs into synth/prior_opt/ (repo left clean).

It re-parses EQUATIONS.txt into shared atoms (reusing the read-only parsing
functions in solve_lab/poly_atoms.py), then:

  1. classifies atoms: unit pins (x-1), booleans (x^2-x), pure copies (x_a-x_b),
     big-constant "load" atoms;
  2. runs integer constraint propagation to fixpoint (reusing the read-only
     solver in solve_lab/propagate.py) WITHOUT touching the 256 core bits, to
     count how many of the 38,748 unknowns are DETERMINED before any search;
  3. builds copy-classes (union-find) over pure copy atoms and locates the class
     of x_26064, checking whether any of the 256 core selector bits are wire-tied;
  4. checks whether any of the 256 core bit variables get forced by propagation
     (they must not -- they are the genuine free variables);
  5. cross-validates SOUNDNESS: every forced value must agree with the best
     known 39,026/39,033 witness (a concrete solution of the non-core part), so
     fixing the forced vars removes no genuine solution.

Outputs: forced_report.json  (+ printed summary).
"""
import json, os, sys, time
from collections import defaultdict, deque

HERE = os.path.dirname(os.path.abspath(__file__))
LAB  = os.path.abspath(os.path.join(HERE, '..', '..', '..'))   # solve_lab
ROOT = os.path.dirname(LAB)
sys.path.insert(0, LAB)

import poly_atoms as PA          # read-only import (guarded by __main__)
import propagate as PR           # read-only import (guarded by __main__)
import ast

NVARS = 38748
EQ_PATH = os.path.join(ROOT, 'EQUATIONS.txt')

# ---------------------------------------------------------------------------
# 1. Re-parse EQUATIONS.txt -> atoms (canon poly dict), in memory.
# ---------------------------------------------------------------------------
def build_atoms():
    atom_eqs = defaultdict(list)
    n = 0
    with open(EQ_PATH) as f:
        for i, line in enumerate(f):
            line = line.strip()
            if not line:
                continue
            n += 1
            lhs = line.rsplit('=', 1)[0]
            core = PA.strip_outer(ast.parse(lhs, mode='eval').body)
            for term in PA.flatten_add(core):
                atom = PA.strip_coef(term)
                poly = PA.expand(atom)
                key = PA.canon(poly)
                if not key:
                    continue
                atom_eqs[key].append(i)
    return atom_eqs, n

def key_to_poly(key):
    return {tuple(m): c for m, c in key}

# ---------------------------------------------------------------------------
def main():
    t0 = time.time()
    atom_eqs, n_eqs = build_atoms()
    keys = list(atom_eqs.keys())
    atoms = [key_to_poly(k) for k in keys]
    print(f"[1] parsed {n_eqs} equations -> {len(atoms)} distinct atoms in {time.time()-t0:.1f}s")

    # classify -------------------------------------------------------------
    unit_pins = []      # x = 1  (atom x-1 : {():-1,(v,):1})
    zero_pins = []      # x = 0  (atom just {(v,):c})
    booleans  = []      # x^2 - x
    copies    = []      # x_a - x_b
    bigload   = 0
    for k in keys:
        d = dict(k)
        allv = set()
        for m, _ in k:
            allv.update(m)
        # big constant load atom
        if any(abs(c) >= 10**20 for _, c in k):
            bigload += 1
        if len(allv) == 1:
            v = next(iter(allv))
            # x - 1  or  1 - x
            if d.get((), 0) != 0 and d.get((v,), 0) != 0 and (v, v) not in d:
                # linear single-var: c1*v + c0 = 0
                c1 = d[(v,)]; c0 = d[()]
                if c0 % c1 == 0 and (-c0) // c1 == 1:
                    unit_pins.append(v)
            if d.get((v, v)) == 1 and d.get((v,)) == -1 and () not in d:
                booleans.append(v)
            if len(k) == 1 and list(k)[0][0] == (v,):
                zero_pins.append(v)
        elif len(allv) == 2 and len(k) == 2:
            (m0, c0), (m1, c1) = sorted(k)
            if len(m0) == 1 and len(m1) == 1 and c0 == -c1:  # x_a - x_b
                copies.append((m0[0], m1[0]))
    print(f"[1] unit pins (x=1): {len(unit_pins)}   booleans x^2-x: {len(booleans)}"
          f"   pure copies x_a-x_b: {len(copies)}   big-const load atoms: {bigload}")

    # ---------------------------------------------------------------------
    # 2. Integer propagation to fixpoint (no core bits touched).
    # ---------------------------------------------------------------------
    var_atoms = defaultdict(list)
    for ai, poly in enumerate(atoms):
        for v in PR.atom_vars(poly):
            var_atoms[v].append(ai)

    val = [None] * NVARS
    contradictions = []
    wl = deque(range(len(atoms)))
    inwl = [True] * len(atoms)

    def assign(v, x):
        if val[v] is not None:
            if val[v] != x:
                contradictions.append((v, val[v], x))
            return
        val[v] = x
        for ai in var_atoms[v]:
            if not inwl[ai]:
                inwl[ai] = True
                wl.append(ai)

    while wl:
        ai = wl.popleft(); inwl[ai] = False
        poly = PR.substitute(atoms[ai], val)
        uv = PR.atom_vars(poly)
        if len(uv) == 1:
            tag, res = PR.solve_single(poly)
            if tag == 'val':
                assign(res[0], res[1])
            elif tag == 'contradiction':
                contradictions.append(('atom', ai, None))
    forced = [v for v in range(NVARS) if val[v] is not None]
    print(f"[2] propagation fixpoint: {len(forced)} of {NVARS} variables DETERMINED "
          f"in {time.time()-t0:.1f}s   contradictions={len(contradictions)}")
    n_one = sum(1 for v in forced if val[v] == 1)
    n_zero = sum(1 for v in forced if val[v] == 0)
    n_other = len(forced) - n_one - n_zero
    print(f"    forced values: =1:{n_one}  =0:{n_zero}  other:{n_other}")

    # ---------------------------------------------------------------------
    # 3. Copy classes (union-find) over pure copies.
    # ---------------------------------------------------------------------
    parent = {}
    def find(a):
        parent.setdefault(a, a)
        while parent[a] != a:
            parent[a] = parent[parent[a]]; a = parent[a]
        return a
    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb: parent[ra] = rb
    for a, b in copies:
        union(a, b)
    classes = defaultdict(list)
    for v in list(parent):
        classes[find(v)].append(v)
    # class of x_26064
    cls_p = sorted(classes[find(26064)]) if 26064 in parent else []
    print(f"[3] copy classes: {len(classes)} classes over {len(parent)} vars; "
          f"class(x_26064) size = {len(cls_p)}")

    # ---------------------------------------------------------------------
    # 4. Core bits: are any of the 256 forced or copy-linked to each other?
    # ---------------------------------------------------------------------
    core = json.load(open(os.path.join(HERE, '..', '..', 'core.json')))
    bit_vars = list(core['bit_vars'])
    core_set = set(bit_vars)
    forced_core = [b for b in bit_vars if val[b] is not None]
    # copy-links among core bits
    core_class = defaultdict(list)
    for b in bit_vars:
        if b in parent:
            core_class[find(b)].append(b)
    linked_pairs = [grp for grp in core_class.values() if len(grp) > 1]
    print(f"[4] of 256 core bits: forced by propagation = {len(forced_core)}; "
          f"copy-linked groups among core bits = {len(linked_pairs)}")

    # ---------------------------------------------------------------------
    # 5. Soundness cross-check vs best 39,026 witness.
    # ---------------------------------------------------------------------
    wit = json.load(open(os.path.join(LAB, 'best', 'new_instance_partial_39026.json')))
    witv = {}
    for kk, vv in wit.items():
        idx = int(kk[2:]) if str(kk).startswith('x_') else int(kk)
        witv[idx] = int(vv)
    # witness omits zeros; forced==0 must be absent-or-zero in witness
    disagree = []
    for v in forced:
        wv = witv.get(v, 0)
        if wv != val[v]:
            disagree.append((v, val[v], wv))
    print(f"[5] soundness vs 39,026 witness: forced vars checked={len(forced)}  "
          f"disagreements={len(disagree)}")
    if disagree[:10]:
        print("    sample disagreements (var, forced, witness):", disagree[:10])

    # ---------------------------------------------------------------------
    report = {
        'n_eqs': n_eqs,
        'n_vars_total': NVARS,
        'n_atoms': len(atoms),
        'unit_pins': len(unit_pins),
        'booleans': len(booleans),
        'pure_copies': len(copies),
        'bigload_atoms': bigload,
        'forced_total': len(forced),
        'forced_eq_1': n_one,
        'forced_eq_0': n_zero,
        'forced_other': n_other,
        'contradictions': len(contradictions),
        'copy_classes': len(classes),
        'class_x26064_size': len(cls_p),
        'core_bits': len(bit_vars),
        'core_bits_forced': len(forced_core),
        'core_bits_copylinked_groups': len(linked_pairs),
        'core_bits_copylinked_detail': linked_pairs,
        'soundness_disagreements': len(disagree),
        'soundness_sample': disagree[:20],
        'free_after_shrink': NVARS - len(forced),
    }
    json.dump(report, open(os.path.join(HERE, 'forced_report.json'), 'w'), indent=1)
    # also dump the forced assignment for an independent checker run
    fa = {f"x_{v}": val[v] for v in forced}
    json.dump(fa, open(os.path.join(HERE, 'forced_assignment.json'), 'w'))
    print(f"[6] wrote forced_report.json ({len(forced)} forced) and forced_assignment.json")
    print(f"    total time {time.time()-t0:.1f}s")

if __name__ == '__main__':
    main()
