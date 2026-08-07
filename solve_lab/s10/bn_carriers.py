"""bn_carriers: operational classification of boolean variables as carriers.

'clean carrier' = flipping it to a non-boolean value makes ONLY its own boolean
atom nonzero (zero downstream damage), so its cost is exactly |E(a)|.
"""
import os, sys, json, collections
HERE='/home/user/integer_solver/solve_lab/s10'
sys.path.insert(0, HERE)
import bn_lib as B
import lib as L, tools as T, ad

bools=B.bools_map()
rows=[json.loads(l) for l in open(os.path.join(HERE,'bn_sweep.jsonl'))]
clean=[r for r in rows if r['nzatoms']==8]           # 7 baseline + the boolean atom
print(f'free boolean vars swept: {len(rows)}')
print(f'CLEAN carriers (only own atom becomes nonzero): {len(clean)}')
h=collections.Counter(39026-r['score'] for r in clean)
print('cost histogram for clean carriers:',sorted(h.items()))
best=min(clean,key=lambda r:39026-r['score'])
print('cheapest clean carrier:',best,'|E(a)|=',len(L.atom2eq[best['atom']]))
# near-clean
h2=collections.Counter(r['nzatoms'] for r in rows)
print('nzatoms histogram (8 = clean):',sorted(h2.items()))
h3=collections.Counter(39026-r['score'] for r in rows)
print('overall cost histogram:',sorted(h3.items())[:14])
print('MIN COST over all 1156 single free-boolean carriers:',min(39026-r['score'] for r in rows))
# syntactic count of "other atoms" for clean carriers
h4=collections.Counter(len(L.var_atoms[r['var']])-1 for r in clean)
print('syntactic other-atom count for clean carriers:',sorted(h4.items()))
json.dump({'clean':[r['var'] for r in clean],
           'min_cost':min(39026-r['score'] for r in rows)},
          open(os.path.join(HERE,'bn_carriers.json'),'w'))
print('saved bn_carriers.json')
