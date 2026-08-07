import os, sys, collections, json
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L
for a in [2202,16897,21113,38521,39166,40066,40932,29539,37887,1459,8261,40005,40121,29090,22231]:
    g=L.atom_out.get(a)
    print(f'a{a} neq={len(L.atom2eq[a])} gate_out={g} nvars={len(L.avars[a])}')
    print(f'   {L.atom_src[a][:400]}')
    print()
