#!/usr/bin/env python3
"""Characterize the ACTIVE-wire-product set: for each wire member, in exact Z at the best partial,
is it used against a NONZERO partner (product) or as a standalone term in a non-identity atom?
The old instance had ZERO active. Identify the minimal active set for the new instance."""
import json, re
from collections import defaultdict
from agentE_common import build_wire, load_best, load_gates, p, NVARS, CORE, NONCORE13

wire, find2, A = build_wire()
best = load_best()
def V(v): return best.get(v, 0)

# classify each atom that touches a wire member
# roles: 'identity' (2-term alias defining wire), 'product' (wire in a degree-2 monomial),
#        'linear' (wire as standalone linear term in a >2-term atom), 'loadcheck' (x_26064 single)
wire_set = set(wire)
active_products = []   # (wire_member, partner, partner_val, atom_idx)
active_linear = []     # (wire_member, atom_idx, natoms_size)
identity_atoms = 0

for ai, poly in enumerate(A):
    vs = atom_vars_set = set()
    for m in poly:
        for v in m: atom_vars_set.add(v)
    wmembers = atom_vars_set & wire_set
    if not wmembers: continue
    # is it a pure 2-term identity?
    keys = list(poly.keys())
    nonconst = [m for m in keys if m != ()]
    is_ident = (len(vs:=atom_vars_set) == 2 and poly.get((), 0) == 0 and
                all(len(m) == 1 for m in nonconst))
    if is_ident:
        # check |c|==|c|
        vv = sorted(atom_vars_set)
        c1 = poly.get((vv[0],), 0); c2 = poly.get((vv[1],), 0)
        if abs(c1) == abs(c2) and c1 and c2:
            identity_atoms += 1
            continue
    # not a defining identity: examine each monomial for wire involvement
    for m, c in poly.items():
        wm = [v for v in m if v in wire_set]
        if not wm: continue
        if len(m) == 2:
            # product monomial: wire * partner  (or wire*wire)
            a, b = m
            if a in wire_set and b in wire_set:
                # wire*wire (e.g. x_26064^2 checks) -- partner is wire
                active_products.append(('WW', a, b, V(a), V(b), ai))
            else:
                partner = b if a in wire_set else a
                wmem = a if a in wire_set else b
                active_products.append(('WP', wmem, partner, V(partner), None, ai))
        elif len(m) == 1:
            active_linear.append((m[0], ai, len(nonconst)))
        # len>2 shouldn't happen (degree<=2 mostly, deg4 are squares)

# summarize product partners with NONZERO value (these are the truly active ones)
print(f"=== identity atoms defining wire: {identity_atoms} ===")
print(f"\n=== ACTIVE PRODUCTS (wire in a degree-2 monomial) ===")
nz_partners = defaultdict(set)   # partner -> set of wire members
ww = []
for rec in active_products:
    if rec[0] == 'WW':
        ww.append(rec)
    else:
        _, wmem, partner, pval, _, ai = rec
        if pval != 0:
            nz_partners[partner].add(wmem)
print(f"wire*wire monomials: {len(ww)} (atoms: {sorted(set(r[5] for r in ww))[:10]})")
print(f"nonzero-partner products: {len(nz_partners)} distinct partners")
for partner, wmems in sorted(nz_partners.items(), key=lambda x: -abs(V(x[0]))):
    print(f"  partner x_{partner} (val bits={V(partner).bit_length()}, sign={'-' if V(partner)<0 else '+'}) <- wire members {sorted(wmems)}")

# zero-partner products (inert at best but would activate if partner moves)
zero_partners = set()
for rec in active_products:
    if rec[0] == 'WP' and rec[3] == 0:
        zero_partners.add(rec[2])
print(f"\nzero-valued partners (inert now): {len(zero_partners)}")

print(f"\n=== ACTIVE LINEAR (wire member standalone in a multi-term atom) ===")
lin_by_member = defaultdict(list)
for wmem, ai, sz in active_linear:
    lin_by_member[wmem].append((ai, sz))
print(f"{len(lin_by_member)} wire members appear as standalone linear terms in non-identity atoms")
for wmem, lst in sorted(lin_by_member.items()):
    print(f"  x_{wmem}: in {len(lst)} atoms, e.g. {lst[:3]}")
