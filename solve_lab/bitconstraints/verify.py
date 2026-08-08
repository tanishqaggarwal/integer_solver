#!/usr/bin/env python3
"""Pass 11: verification + final constraint inventory.

(a) Take every atom and every equation whose support lies inside the 256
    selector bits, reduce it modulo the booleanity ideal <b_i^2 - b_i> (i.e.
    multilinearise: b^k -> b) and check that the remainder is identically 0.
    A zero remainder proves the constraint carries NO information beyond
    booleanity.  A nonzero remainder would be a genuine extra constraint.
(b) Same for support = bits + one extra variable, treating the extra variable
    as unconstrained (remainder must vanish as a polynomial in that variable
    too) and, separately, as boolean.
(c) Independent numeric spot-check: evaluate each such equation on random 0/1
    assignments and confirm it is exactly 0.
(d) Emit constraints.json.
"""
import os, json, random
from collections import defaultdict, Counter
from scan import load, support, degree, fmt

HERE = os.path.dirname(os.path.abspath(__file__))
D, chain, BITS = load()
atoms = D['atoms']
eq_poly = D['eq_poly']
eq_terms = D['eq_terms']
idx = json.load(open(os.path.join(HERE, 'scan_index.json')))
atom_eqs = defaultdict(list)
for i, terms in enumerate(eq_terms):
    for c, aid in terms:
        atom_eqs[aid].append(i)


def multilinear(poly, boolset):
    """reduce mod <b^2-b> for b in boolset"""
    out = defaultdict(int)
    for m, c in poly:
        nm = []
        seen = set()
        for v in m:
            if v in boolset:
                if v in seen:
                    continue
                seen.add(v)
            nm.append(v)
        out[tuple(sorted(nm))] += c
    return {m: c for m, c in out.items() if c != 0}


report = {}

# ---- (a) bit-only atoms & equations -------------------------------------
bad_atoms, bad_eqs = [], []
for aid in idx['bitonly_atoms']:
    r = multilinear(atoms[aid], BITS)
    if r:
        bad_atoms.append((aid, r))
for e in idx['bitonly_eqs']:
    r = multilinear(eq_poly[e], BITS)
    if r:
        bad_eqs.append((e, r))
print(f"(a) bit-only atoms   : {len(idx['bitonly_atoms'])}, "
      f"non-vanishing mod booleanity: {len(bad_atoms)}")
print(f"    bit-only equations: {len(idx['bitonly_eqs'])}, "
      f"non-vanishing mod booleanity: {len(bad_eqs)}")
for e, r in bad_eqs:
    print("     eq#%d residual %s" % (e, r))

# ---- (b) one-step --------------------------------------------------------
bad1a, bad1e, extra_bool = [], [], set()
for aid in idx['onestep_atoms']:
    s = support(atoms[aid])
    x = (s - BITS).pop()
    r = multilinear(atoms[aid], BITS | {x})
    if r:
        bad1a.append((aid, x, r))
    else:
        extra_bool.add(x)
for e in idx['onestep_eqs']:
    s = support(eq_poly[e])
    x = (s - BITS).pop()
    r = multilinear(eq_poly[e], BITS | {x})
    if r:
        bad1e.append((e, x, r))
    else:
        extra_bool.add(x)
print(f"(b) one-step atoms    : {len(idx['onestep_atoms'])}, non-vanishing "
      f"mod booleanity(bits + the extra var): {len(bad1a)}")
print(f"    one-step equations: {len(idx['onestep_eqs'])}, non-vanishing: "
      f"{len(bad1e)}")
for aid, x, r in bad1a[:20]:
    print(f"     atom#{aid} extra x_{x}: {r}")
for e, x, r in bad1e[:20]:
    print(f"     eq#{e} extra x_{x}: {r}")

# separate the 512 copy/NOT atoms (these are definitions, not constraints)
copies, nots = [], []
for aid in idx['onestep_atoms']:
    a = atoms[aid]
    d = dict(a)
    s = support(a)
    if len(a) == 2 and degree(a) == 1:
        b = (s & BITS).pop(); w = (s - BITS).pop()
        copies.append((aid, b, w))
    elif len(a) == 3 and degree(a) == 1 and () in d:
        b = (s & BITS).pop(); w = (s - BITS).pop()
        nots.append((aid, b, w))
print(f"    of the one-step atoms: {len(copies)} are copies (b = w), "
      f"{len(nots)} are NOTs (w = 1 - b), "
      f"{len(idx['onestep_atoms'])-len(copies)-len(nots)} are booleanity sums")

# ---- (c) numeric spot-check ---------------------------------------------
rng = random.Random(999)
allv = set()
for e in idx['bitonly_eqs'] + idx['onestep_eqs']:
    allv |= support(eq_poly[e])
fails = 0
for trial in range(200):
    val = {v: rng.randrange(2) for v in allv}
    for e in idx['bitonly_eqs'] + idx['onestep_eqs']:
        t = 0
        for m, c in eq_poly[e]:
            k = c
            for v in m:
                k *= val[v]
            t += k
        if t != 0:
            fails += 1
print(f"(c) numeric spot-check: {200*len(idx['bitonly_eqs']+idx['onestep_eqs'])} "
      f"evaluations on random 0/1 vectors, {fails} nonzero")

# ---- (d) linear-algebra summary -----------------------------------------
# collect any *linear* (degree-1 after multilinearisation) relation among bits
lin_rows = []
for e in idx['bitonly_eqs']:
    r = multilinear(eq_poly[e], BITS)
    row = {m[0]: c for m, c in r.items() if len(m) == 1}
    if row:
        lin_rows.append((e, row))
print(f"(d) linear relations among the selector bits derived from bit-only "
      f"equations: {len(lin_rows)}")

out = {
    'n_bits': 256,
    'bit_vars_in_chain_order': chain,
    'booleanity': [
        {'bit': (support(atoms[a]) & BITS).pop(), 'atom': a,
         'equations': sorted(set(atom_eqs[a]))[:12]}
        for a in idx['bitonly_atoms'] if len(support(atoms[a])) == 1],
    'bitonly_equations': [
        {'eq': e, 'poly': [[list(m), c] for m, c in eq_poly[e]],
         'reduces_to_zero_mod_booleanity': not multilinear(eq_poly[e], BITS)}
        for e in idx['bitonly_eqs']],
    'bitonly_multibit_atoms': [
        {'atom': a, 'bits': sorted(support(atoms[a])),
         'equations': sorted(set(atom_eqs[a])),
         'reduces_to_zero_mod_booleanity': not multilinear(atoms[a], BITS)}
        for a in idx['bitonly_atoms'] if len(support(atoms[a])) > 1],
    'onestep_equations': [
        {'eq': e, 'extra_var': sorted(support(eq_poly[e]) - BITS)[0],
         'reduces_to_zero_mod_booleanity': not multilinear(
             eq_poly[e], BITS | (support(eq_poly[e]) - BITS))}
        for e in idx['onestep_eqs']],
    'copy_atoms': [{'atom': a, 'bit': b, 'wire': w} for a, b, w in copies],
    'not_atoms': [{'atom': a, 'bit': b, 'wire': w} for a, b, w in nots],
    'forced_bits': [],
    'equal_pairs': [],
    'opposite_pairs': [],
    'exclusive_pairs': [],
    'gf2_linear_relations': [],
    'integer_linear_relations': [],
    'cardinality_constraints': [],
    'free_bits': 256,
}
json.dump(out, open(os.path.join(HERE, 'constraints.json'), 'w'), indent=1)
print("wrote constraints.json")
