import os, sys, collections
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad
P = 2**256-2**32-977
for a in (37887, 7930, 29539, 41512, 40826, 21617):
    try:
        print(f'=== a{a}  eqs={sorted(L.atom2eq.get(a,{}).items())}  out={L.atom_out.get(a)}')
        print('SRC:', L.atom_src[a])
        print()
    except Exception as ex: print(a, ex)
