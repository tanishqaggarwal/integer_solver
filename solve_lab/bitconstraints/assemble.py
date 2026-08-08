#!/usr/bin/env python3
"""Pass 16: assemble the machine-readable constraint inventory."""
import os, json
from collections import defaultdict
from scan import load, support, degree, fmt

HERE = os.path.dirname(os.path.abspath(__file__))
D, chain, BITS = load()
atoms = D['atoms']
eq = D['eq_poly']
idx = json.load(open(os.path.join(HERE, 'scan_index.json')))
ortree = json.load(open(os.path.join(HERE, 'ortree.json')))
bp = json.load(open(os.path.join(HERE, 'booleanity_proof.json')))
sweep = json.load(open(os.path.join(HERE, 'sweep.json')))
atom_eqs = defaultdict(list)
for i, terms in enumerate(D['eq_terms']):
    for c, aid in terms:
        atom_eqs[aid].append(i)


def poly_json(p):
    return [[list(m), c] for m, c in p]


bool_atoms = {}
for aid in idx['bitonly_atoms']:
    s = support(atoms[aid])
    if len(s) == 1:
        bool_atoms[next(iter(s))] = aid

orc = ortree['or_constraints'][0]

out = {
  "instance": {
    "equations_file": "EQUATIONS.txt",
    "n_equations": len(eq),
    "n_distinct_atoms": len(atoms),
    "n_selector_bits": 256,
    "selector_bit_vars_chain_order": chain,
  },
  "summary": {
    "forced_bits": [],
    "equal_pairs": [],
    "opposite_pairs": [],
    "mutually_exclusive_pairs": [],
    "gf2_linear_relations": [],
    "integer_linear_relations": [],
    "rank_of_linear_constraint_space_on_bits": 0,
    "free_bits": 256,
    "search_space_before": "2**256",
    "search_space_after": "2**256 - 1   (only the all-zero vector is excluded)",
  },
  "booleanity": {
    "description": "b*(b-1) = 0 for each of the 256 selectors",
    "atoms": [{"bit": b, "atom": a,
               "equations": sorted(set(atom_eqs[a]))}
              for b, a in sorted(bool_atoms.items())],
    "rigorous_sign_certificates": {
      "method": "for integer x, x - x^2 <= 0 with equality iff x in {0,1}; an "
                "equation that is a sign-uniform Z-combination of such terms "
                "forces every term to 0",
      "n_pure_booleanity_combination_equations": bp['n_bool_combination_equations'],
      "vars_covered": bp['covered_vars'],
      "selectors_proven_boolean": bp['proven_boolean_selectors'],
      "certificate_equation_per_var": bp['certificate_equation_per_var'],
    },
  },
  "cardinality": [{
    "constraint": "OR(b_0, ..., b_255) = 1   (at least one selector set; k != 0)",
    "n_selectors_in_the_OR": orc['n_selectors'],
    "selectors": orc['selectors'],
    "root_wire_class": orc['wire_class'],
    "root_pins": orc['pins'],
    "or_gadget_atoms": orc['gadget_atoms'],
    "gadget_shape": "o = u OR v encoded by three atoms: (s - u - v), (u*v - p), (s - p - o)",
    "n_or_gadgets": ortree['n_or_gadgets'],
    "note": "implied by the accumulator condition (T != O), and separately "
            "encoded in the circuit; it is the ONLY pure-bit constraint the "
            "system adds on top of booleanity",
  }],
  "bitonly_equations": [
    {"eq": e,
     "support": sorted(support(eq[e])),
     "poly": poly_json(eq[e]),
     "is_in_the_booleanity_ideal": True,
     "linear_part_after_multilinearisation": []}
    for e in idx['bitonly_eqs']],
  "bitonly_multibit_atoms": [
    {"atom": a, "bits": sorted(support(atoms[a])),
     "equations": sorted(set(atom_eqs[a])),
     "poly": poly_json(atoms[a]),
     "is_in_the_booleanity_ideal": True}
    for a in idx['bitonly_atoms'] if len(support(atoms[a])) > 1],
  "onestep_equations": [
    {"eq": e, "extra_var": sorted(support(eq[e]) - BITS)[0],
     "bits": sorted(support(eq[e]) & BITS),
     "is_in_the_booleanity_ideal_of_bits_plus_extra": True}
    for e in idx['onestep_eqs']],
  "onestep_atoms_are_definitions": {
    "copy_atoms": "256 atoms  b - w   (alias wire)",
    "not_atoms": "256 atoms  1 - b - w  (complement wire)",
    "booleanity_sums": 6,
    "conclusion": "no new constraint on the bits; they only name 512 derived wires",
  },
  "pin_gates": {
    "count": 512,
    "per_bit": 2,
    "shape": "b*w - C*b - s*v = 0,  i.e.  b*(w - C) = s*v",
    "collisions": 0,
    "note": "the only other way a selector enters the system; no two selectors "
            "gate the same target with different constants, so no bit is forced "
            "and no pair is excluded",
  },
  "local_consistency_sweep": {
    "n_patterns_tested": sweep['n_cases'],
    "patterns": "all-zeros, all-ones, all 256 weight-1, 40 weight-2, 40 random",
    "n_primitive_atoms_used": sweep['n_prim_atoms'],
    "patterns_with_a_propagation_contradiction": sweep['contradicting_patterns'],
  },
}
json.dump(out, open(os.path.join(HERE, 'constraints.json'), 'w'), indent=1)
print("wrote constraints.json",
      os.path.getsize(os.path.join(HERE, 'constraints.json')))
