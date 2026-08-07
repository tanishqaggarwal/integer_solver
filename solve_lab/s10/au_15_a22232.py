import os, sys, collections, itertools
from fractions import Fraction
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad
P = 2**256-2**32-977
SEVEN=[22229,22230,35758,35759,35760,35761,35762]
E12=[2554,6816,8124,9123,9421,12231,12270,12350,14584,18673,22044,29125]
QATOMS=[(1,22231),(6,22232),(15,22233),(-21,22234),(-13,22235),(25,19087),(1,19088),
        (25,19089),(28,19090),(1,19091),(-4,19092),(23,10935),(-5,10936),(20,10937),
        (-27,10938),(35,10939),(17,10940),(-14,10941)]
v=L.load(os.path.join(LAB,'best','new_instance_partial_39026.json'))
av=L.all_atom_values(v)
print('outside-12 equations for each Q atom:')
for c,a in QATOMS:
    eqs=set(L.atom2eq[a]); out=sorted(eqs-set(E12))
    print(f'  a{a:<6} coeffQ={c:<5} n_eq={len(eqs)} in12={len(eqs&set(E12))} out={len(out)} {out}')
print()
# For a22232: what is its single outside equation, and can it be satisfied?
for a in (22232,):
    out=sorted(set(L.atom2eq[a])-set(E12))
    for e in out:
        m,sq,co=L.eq_atoms[e]
        print(f'a{a} outside eq {e}: mult={m} sq={sq} natoms={len(co)}')
        for aa,cc in sorted(co.items()):
            print(f'    a{aa:<6} c={cc:<5} val_nonzero={av[aa]!=0} n_eqs={len(L.atom2eq[aa])} out_of_12={len(set(L.atom2eq[aa])-set(E12))}')
