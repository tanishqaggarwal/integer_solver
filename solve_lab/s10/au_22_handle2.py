import os, sys, collections
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad
P = 2**256-2**32-977
v=L.load(os.path.join(LAB,'best','new_instance_partial_39026.json'))
for u in (38744, 3387, 22972, 5081, 11627, 17286, 34817, 13632):
    d=L.definer.get(u)
    print(f'x_{u}: free={d is None} definer=a{d} val={str(v[u])[:50]} ==p:{v[u]==P} bits={v[u].bit_length()} n_atoms={len(L.var_atoms[u])}')
    if d is not None: print('     def src:', L.atom_src[d][:150])
print()
C6418=33472904810391811973223207617762334363023286939839396241234196646906030803538671321618319
C12553=42775533402728869434716629464193396056515231264222641773817154079369026410240838606908039
print('x_6418  =', v[6418]); print('C6418   =', C6418); print('diff =', v[6418]-C6418)
print('x_12553 =', v[12553]); print('C12553  =', C12553); print('diff =', v[12553]-C12553)
print('C6418 mod p =', C6418%P, ' x_6418 mod p =', v[6418]%P)
print('C6418 // p =', C6418//P, ' 15804267*p =', 15804267*P)
print('diff/15804267 =', (v[6418]-C6418)/15804267 if (v[6418]-C6418)%15804267==0 else 'not divisible')
