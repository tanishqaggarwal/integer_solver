import os, sys, collections, json
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad
P = 2**256-2**32-977
SEVEN = [22229, 22230, 35758, 35759, 35760, 35761, 35762]
v = L.load(os.path.join(LAB,'best','new_instance_partial_39026.json'))
av = L.all_atom_values(v)
for a in SEVEN + [22231, 22232, 22233, 22234, 22235, 35757, 35756, 35755, 37887, 7930, 29539]:
    print(f'a{a}: out={L.atom_out.get(a)}  src={L.atom_src[a][:220]}')
    print(f'     value={av[a]}   value mod p = {av[a]%P}')
    print(f'     vars={sorted(L.avars[a])}')
    print()
