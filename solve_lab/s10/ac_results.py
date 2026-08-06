"""S12: consolidate every measurement of the activation study into one durable
JSON (ac_results.json).  Pure bookkeeping -- no computation."""
import os, json
HERE='/home/user/integer_solver/solve_lab/s10'
def g(p):
    try: return json.load(open(os.path.join(HERE,p)))
    except Exception: return None
R = {
 'frames': {
   'f39009': {'file':'s10/mod9118_0.json','score':39009,
     'nonzero_atoms':[21617,29539,37662,40826],'eqs_touched':24,
     'cluster_cone_vars':1401,'zero_free_inputs_in_cone':189,
     'grad_support_of_cluster':82,'ceiling_if_21617_29539_zeroed':39031},
   'f39026': {'file':'best/new_instance_partial_39026.json','score':39026,
     'nonzero_atoms':[22229,22230,35758,35759,35760,35761,35762],'eqs_touched':12,
     'cluster_cone_vars':39,'zero_free_inputs_in_cone':1,
     'grad_support_of_cluster':3,
     'note':'the record frame is structurally immune to activation'}},
 'singles_39009': {'candidates':189,'values_tested':[1,-1,'generic'],
   'zero_equation_cost':42,'grow_support':76,'max_knobs':2,
   'best':{'x':24365,'val':1,'knobs':2,'eqs_lost':19,'checks_broken':5,'score':38990}},
 'pairs_39009_val1': {'pairs':17766,'grew_support':11438,'max_knobs':2,
   'genuine_second_order':0,'zero_cost_growth':0,
   'best':{'x':[1081,24365],'knobs':2,'eqs_lost':19,'score':38990}},
 'pairs_39009_generic': {'pairs':17766,'grew_support':11438,'max_knobs':4,
   'genuine_second_order':1846,
   'best':{'x':[11368,17406],'knobs':4,'eqs_lost':71,'score':38938}},
 'second_py_correction':
   'second.py reported 6/6 pairs growing the support; all six partners '
   '(x_12054,x_16586,x_17406,x_28713,x_27393,x_11368) grow it as SINGLES, so the '
   'reported gain belongs to z alone, not to the pair.',
 'closure_39009_base': {'rows':1655,'cols':707,'rank':707,'kernel':0,
   'inconsistent_reduced_rows':11,'min_eq_cost_decode':15,'mod_p_ceiling':39018,
   'isd_information_sets':1035,
   'best_violated_set':[21617,37662,41400,41507,41827,42245],
   'drop_all_396_weight1_checks':'still inconsistent',
   'consistent_only_at_weight_ge_12':{'rows':807,'rank':512,'kernel':195,'eq_cost':69}},
 'closure_after_activation': g('ac_closure2.json'),
 'closure_after_generic_activation': g('ac_closure_gen.json'),
 'kernel_answer':'rank == cols in every targeted configuration; the only kernel '
   'ever seen is 2/976 at full 189-input generic activation, with 361 inconsistent '
   'rows and a 962-equation decode cost.',
 'rowcost': {'measured':152,'no_new_rows_and_no_lost_cols':73,
   'note':'the new knobs themselves add no closure row, but rebuilding shows rows '
          'still grow because the existing 707 columns move with the point'},
 'frozen_repair': [
   {'x':24365,'val':1,'activated':38990,'repaired':39003,'knobs':2,'alive':True},
   {'x':12054,'val':1,'activated':38989,'repaired':39003,'knobs':2,'alive':True},
   {'x':16586,'val':1,'activated':38988,'repaired':39003,'knobs':2,'alive':True},
   {'x':17406,'val':1,'activated':38988,'repaired':38990,'knobs':2,'alive':True},
   {'x':11368,'val':1,'activated':38980,'repaired':38980,'knobs':2,'alive':True},
   {'x':24365,'val':'gen','activated':38975,'repaired':38975,'knobs':2,'alive':True},
   {'x':12054,'val':'gen','activated':38975,'repaired':38976,'knobs':2,'alive':True},
   {'x':16586,'val':'gen','activated':38975,'repaired':38975,'knobs':2,'alive':True},
   {'x':17406,'val':'gen','activated':38973,'repaired':38973,'knobs':2,'alive':True},
   {'x':11368,'val':'gen','activated':38966,'repaired':38966,'knobs':2,'alive':True}],
 'unfrozen_repair':'the engine repairs the collateral by UNDOING the activation: '
   'engine_ac_a12054_39009.json / a16586 / none are byte-identical to mod9118_0.json '
   '(0 variables differ, score back to 39009).',
 'stacked_frozen': [
   {'k':1,'raw':38990,'repaired':39003,'knobs':2,'closure':[1663,709],'rank':709,
    'kernel':0,'min_eq_cost':19},
   {'k':2,'raw':38978,'repaired':38978,'knobs':2},
   {'k':3,'raw':38975,'repaired':38975,'knobs':2}],
 'sweep_39026': {'single_moves_tested':14541,'free_inputs':7273,
   'best_score':39026,'improving_moves':0,'score_neutral_moves':3878,
   'note':'no single +-1 move over ANY free input beats 39026; all neutral moves '
          'lie outside the 39-variable cluster cone so they can never reach it'},
 'best_score_produced':39009,
 'verdict':'activation never beat 39026; no kernel opened in any targeted '
   'configuration; the record frame has nothing to activate.'}
json.dump(R, open(os.path.join(HERE,'ac_results.json'),'w'), indent=1)
print('wrote ac_results.json')
for k in R: print(' ', k)
