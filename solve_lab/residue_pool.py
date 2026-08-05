#!/usr/bin/env python3
"""Do the 22-side and 233-side draw from the SAME residue pool?
For each control bit, find the HUGE residue it gates via a simple huge-atom
bit*(x_B - HUGE) - s*x_C. Compare the 22-side vs 233-side pools. Also print the
exact LINEAR form of x_8821 (shared denominator) over the 233 bits."""
import json
from collections import defaultdict
from confluent_eval5 import build5, make_forward, boolean_vars
from propagate import load_atoms, atom_vars, NVARS

BITS22 = set([1782,1858,2795,2800,3483,5443,10652,19520,21188,21588,23634,26947,
          27512,29682,30104,30596,30658,30792,33251,37748,37885,38116])

def main():
    A = load_atoms()
    bset = boolean_vars(A)
    control = set(json.load(open('control_bits.json')))

    # each simple huge-atom: 1 product (bit*x_B), 1 huge*bit const term, 1 linear (s*x_C)
    bit_res = defaultdict(list)   # control bit -> list of HUGE residues it gates
    for a, poly in enumerate(A):
        prod2 = [(m, c) for m, c in poly.items() if len(m) == 2]
        huges = [(m, c) for m, c in poly.items() if abs(c) > 10**40]
        if len(prod2) != 1 or not huges: continue
        (m2, c2) = prod2[0]
        for bit in m2:
            if bit in control:
                # HUGE is the coef of the (bit,) linear term OR product; take max-abs huge
                h = max((abs(c) for m, c in huges), default=0)
                bit_res[bit].append(h)
    pool22 = defaultdict(int); pool233 = defaultdict(int)
    for b, rs in bit_res.items():
        for r in rs:
            if b in BITS22: pool22[r] += 1
            else: pool233[r] += 1
    set22 = set(pool22); set233 = set(pool233)
    print(f"22-side bits gating a huge residue: {sum(1 for b in bit_res if b in BITS22)}/22")
    print(f"233-side bits gating a huge residue: {sum(1 for b in bit_res if b not in BITS22)}/233")
    print(f"distinct residues: 22-side pool={len(set22)}, 233-side pool={len(set233)}, SHARED={len(set22&set233)}")
    if set22 & set233:
        print("  shared residues (first 3):", [str(x)[:40]+'...' for x in list(set22&set233)[:3]])
    # residue magnitude distribution
    allres = set22 | set233
    bits_per = sorted(set(len(bin(r))-2 for r in allres))
    print(f"residue bit-lengths present: {bits_per}")

    # exact linear form of x_8821 over the 233 bits (Z, single-flip; linear so exact)
    Aatoms, kind, info, seq0, bestval, ncyc = build5()
    order = json.load(open('eval_order.json'))['order']
    defset = set(v for v in kind if kind[v] != 'const')
    seq = [v for v in order if v in defset and v not in (9770, 3183)]
    seq += [v for v in (9770, 3183) if v in defset]
    seq += [v for v in defset if v not in set(order) and v not in (9770, 3183)]
    solve = make_forward(kind, info, seq, bestval)
    base = solve(list(bestval), [])
    c8 = base[8821]
    bits233 = [b for b in control if b not in BITS22]
    weights = {}
    for b in bits233:
        val = solve(list(bestval), [b])
        d = val[8821] - c8
        if d: weights[b] = d
    print(f"\nx_8821 linear form: base={c8}, {len(weights)} nonzero weights over 233 bits")
    # are the weights themselves huge residues?
    wvals = sorted(set(abs(w) for w in weights.values()))
    print(f"  distinct |weight| values: {len(wvals)}; bit-lengths: {sorted(set(len(bin(w))-2 for w in wvals))}")
    print(f"  sample weights:", [str(weights[b]) for b in list(weights)[:2]])
    # check: do x_8821 weights coincide with the residue pool?
    winset = sum(1 for w in wvals if w in allres)
    print(f"  x_8821 weights that ARE pool residues: {winset}/{len(wvals)}")

    json.dump({'pool22': sorted(str(x) for x in set22),
               'pool233': sorted(str(x) for x in set233),
               'x8821_base': str(c8),
               'x8821_weights': {str(b): str(w) for b, w in weights.items()}},
              open('residue_pool.json', 'w'))
    print("wrote residue_pool.json")

if __name__ == '__main__':
    main()
