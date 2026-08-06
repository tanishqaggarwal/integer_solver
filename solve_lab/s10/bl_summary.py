"""bl_summary: consolidate every measured number into bl_summary.json (durable)."""
import os, sys, json, glob
HERE='/home/user/integer_solver/solve_lab/s10'
S = {}
S['baselines'] = {'witness_frame2_39026': 39026, 'canonical_mod9118_0': 39009}
S['booleans'] = {'total_boolean_vars': 3484, 'free_boolean_inputs': 1156,
                 'free_excluding_MUX_controls': 1154}
S['cones'] = {
  'cone(a21617)': {'vars': 260, 'free_booleans': 10},
  'cone(a29539)': {'vars': 1277, 'free_booleans': 78},
  'cone(a21617,a29539)': {'vars': 1401, 'free_booleans_excl_MUX': 76},
  'cone(a37662)': {'vars': 690, 'free_booleans': 21},
  'cone(a40826)': {'vars': 3217, 'free_booleans': 168},
  'cone(a7930)': {'vars': 74, 'free_booleans': 3, 'which': [2081, 4287, 13195]},
  'frame2_cone(seven_residual_atoms)': {'vars': 39, 'booleans': [2081, 4287],
      'free_booleans_excl_MUX': 0},
  'frame3_cone(eight)': {'vars': 47, 'free_booleans_excl_MUX': 0},
  'cone(all 24 atoms in the 7 failing equations)': {'vars': 151,
      'free_booleans_excl_MUX': [11368, 13195]},
}
S['load_pins'] = {
  'shape': 'b*(x - K) - c*z   (3 monomials)',
  'conditional_constant_loads': 727, 'distinct_atoms': 623,
  'clean_3_monomial_pins': 512, 'distinct_gate_booleans': 256,
  'all_gates_are_free_inputs': True,
  'gates_reaching_C0(x_7068-x_2099)': {2081: [[3576, 6418]], 4287: [[3568, 31861], [36602, 31861]]},
  'gates_reaching_cone(a7930)': {2081: [[3578, 12553]], 4287: [[3570, 14865], [36602, 14865]],
                                 13195: [[3586, 6467]]},
  'gates_reaching_cone(a21617)': [12054, 16586, 17406, 22562, 23751, 24365],
  'gates_reaching_cone(a29539)': [438, 2081, 4287, 5643, 6821, 12054, 13195, 16586,
                                  16827, 17760, 21074, 24365, 27393, 38625],
  'gates_reaching_cone(a37662)': 16, 'gates_reaching_cone(a40826)': 50,
  'pinned_constants_bits': {'a3576(x_2081->x_6418)': 295, 'a3578(x_2081->x_12553)': 295,
     'a3568(x_4287->x_31861)': 296, 'a3570(x_4287->x_14865)': 296,
     'a3584(x_13195->x_18623)': 293, 'a3586(x_13195->x_6467)': 296,
     'a15462(x_11368->x_24559)': 294, 'a15464(x_11368->x_4068)': 294},
}
S['single_flip_scans'] = {
  'frame2_from_39026': {'tested': 1154, 'neutral_at_39026': 900, 'worse': 254, 'better': 0},
  'canonical_from_39009': {'tested': 1154, 'neutral_at_39009': 900, 'worse': 254, 'better': 0},
}
S['mass_flip_test_frame2'] = {
  'rand10': 38975, 'rand50': 38776, 'rand200': 38369, 'all_900_neutral': 39026,
  'random_half_of_neutral(450)': 39026, 'ALL_1154': 36880,
  'seven_residual_atom_values_identical_in_every_case': True}
S['MUX_branch_engine_fixed_points'] = {
  'frame2_from_39026': {'(0,0)': 38984, '(0,1)': 39011, '(1,0)': 39026, '(1,1)': 39003},
  'canonical_from_mod9118_0': {'(0,0)': 38986, '(0,1)': 39009, '(1,0)': 39009, '(1,1)': 38994},
}
S['key4_enumeration'] = {
  'booleans': [2081, 4287, 11368, 13195],
  'note': 'engine fixed point depends only on the MUX pair; x_11368/x_13195 irrelevant',
  'results_frame2': {'(0,0,*,*)': [38984, 39000, 38984, 38984], '(0,1,*,*)': 39011,
                     '(1,0,*,*)': 39026, '(1,1,*,*)': 39003},
  'best': 39026}
S['pair_triple_scans'] = {
  'canonical_pairs_over_76_cluster_cone_booleans': {'combinations': 2850,
      'best_cheap': 38982, 'best_cheap_pair': [16586, 24365], 'base': 39009},
  'priority_scan_canonical': {'cheap_evals': 1232,
      'set': [438,490,5643,5910,6821,12054,12095,13195,14808,16586,16827,17406,
              17760,18022,21074,22562,23751,24365,27393,28005,34974,38625],
      'design': 'singles+pairs over 4 MUX prefixes + 220 triples',
      'best_cheap': 38990, 'engine_shortlist': 14, 'best_engine': 39009, 'base': 39009},
  'branch(1,1)+one_flip_canonical': {'combinations': 1154, 'best_cheap': 38960,
      'prefix_alone': 38960},
  'branch(1,1)+one_flip_frame2': {'combinations': 1154, 'best_cheap': 38983,
      'prefix_alone': 38983},
}
S['conclusion'] = {
  'best_found': 39026, 'improvement_over_baseline': 0, 'bl_best_json_written': False}
if os.path.exists(os.path.join(HERE, 'bl_confirm.jsonl')):
    S['confirm'] = [json.loads(l) for l in open(os.path.join(HERE, 'bl_confirm.jsonl'))]
json.dump(S, open(os.path.join(HERE, 'bl_summary.json'), 'w'), indent=1)
print(json.dumps(S['conclusion']))
print('wrote', os.path.join(HERE, 'bl_summary.json'))
