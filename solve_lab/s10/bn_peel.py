"""bn_peel: which boolean atoms can be nonzero AT ALL?

If an equation e contains exactly one boolean atom a and all its other atoms are
zero, then coeff*val_a = 0 => val_a = 0.  Peel to a fixed point.  What survives
is the maximal possible support of a pure boolean-carrier configuration.

Two regimes:
  strict  - every non-boolean atom is required zero (pure carrier)
  relaxed - allow the 7 already-nonzero atoms of the 39026 state to stay nonzero
"""
import os, sys, json, collections
HERE='/home/user/integer_solver/solve_lab/s10'
sys.path.insert(0, HERE)
import bn_lib as B
import lib as L, tools as T, ad

bools = B.bools_map()
BA = set(bools)
FREEB = set(a for a in BA if bools[a][0] in B.FREESET)

def peel(cand, allow_nz=frozenset()):
    """cand: set of boolean atoms allowed nonzero. allow_nz: non-boolean atoms
    allowed nonzero.  Returns surviving set."""
    cand = set(cand)
    changed = True
    while changed:
        changed = False
        eqs = collections.defaultdict(list)
        for a in cand:
            for e in L.atom2eq[a]: eqs[e].append(a)
        for e, As in eqs.items():
            m, sq, co = L.eq_atoms[e]
            # non-boolean atoms in e that are NOT allowed nonzero -> they are 0
            free_other = any((x not in BA) and (x in allow_nz) for x in co)
            if free_other:
                continue                      # equation has a free absorber
            if len(As) == 1:
                cand.discard(As[0]); changed = True
    return cand

for label, cand in (('all boolean', BA), ('free-var boolean', FREEB)):
    for al, nzset in (('strict', frozenset()),
                      ('relaxed(7 nz atoms free)', frozenset([22229,22230,35758,35759,35760,35761,35762]))):
        s = peel(cand, nzset)
        Es = set()
        for a in s: Es |= set(L.atom2eq[a])
        print(f'{label:20s} {al:26s} survivors={len(s):5d} eqs={len(Es):5d} '
              f'defic={len(Es)-len(s)}')
        if 0 < len(s) <= 40:
            print('    ', sorted(s))
        globals()['S_'+label.split()[0]+'_'+al[:3]] = s

# save the strict all-boolean core
core = peel(BA)
json.dump(sorted(core), open(os.path.join(HERE,'bn_core.json'),'w'))
print('saved bn_core.json', len(core))
# free-only core
fcore = peel(FREEB)
json.dump(sorted(fcore), open(os.path.join(HERE,'bn_fcore.json'),'w'))
print('free core', len(fcore))
