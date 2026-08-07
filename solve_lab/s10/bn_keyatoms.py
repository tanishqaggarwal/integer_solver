"""bn_keyatoms: can the 4 influential bits' boolean atoms be cancelled at all?

These are the ONLY boolean atoms whose variable can reach the 7 failing equations.
For each, list its equations, the other boolean atoms sharing them, and the
per-equation cancellation ratio.  A ratio must be positive AND consistent across
all of the atom's equations for zero-cost activation.
"""
import os, sys, json, collections
from fractions import Fraction
HERE='/home/user/integer_solver/solve_lab/s10'
sys.path.insert(0, HERE)
import bn_lib as B
import lib as L, tools as T, ad

bools=B.bools_map()
inf=json.load(open(os.path.join(HERE,'bn_infl.json')))
KEY=sorted(set(inf['bool_anc2']))
peel=set(json.load(open(os.path.join(HERE,'bn_cone.json')))['blk29_peel'])
byvar={u:a for a,(u,c) in bools.items()}
print('influential bits and their boolean atoms:')
for u in KEY:
    a=byvar[u]
    print(f'  x_{u} -> a{a}  c={bools[a][1]}  |E|={len(L.atom2eq[a])} '
          f'in_maximal_support={a in peel}  free={u in B.FREESET}')
print()
for u in KEY:
    a=byvar[u]
    print(f'=== a{a} (x_{u}) ===')
    for e,co in sorted(L.atom2eq[a].items()):
        m,sq,cc=L.eq_atoms[e]
        others=[(b,cc[b]) for b in cc if b in bools and b!=a]
        oth_live=[(b,c) for b,c in others if b in peel]
        # ratio needed for a partner b to cancel a alone
        rr=[]
        for b,c in oth_live:
            r=Fraction(-cc[a]*bools[a][1], c*bools[b][1])
            rr.append((b, str(r), bools[b][0] in B.FREESET))
        print(f'  eq{e}: natoms={len(cc)} coeff_a={cc[a]} boolean_partners={len(others)} '
              f'live={len(oth_live)}')
        if rr: print(f'     t_b/t_a ratios needed: {rr[:6]}')
        else:  print(f'     NO live boolean partner -> this equation CANNOT be cancelled')
    print()
