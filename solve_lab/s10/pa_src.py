import os, sys, collections, json
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad
S=[22229, 22230, 35758, 35759, 35760, 35761, 35762]
for a in S:
    print(f'a{a}: {L.atom_src[a][:220]}')
print()
E=sorted(set().union(*[set(L.atom2eq[a]) for a in S]))
print('E =',E)
for e in E:
    m,sq,co = L.eq_atoms[e]
    print(f'eq{e}: mult={m} sq={sq} natoms={len(co)}  co_on_S={{ {", ".join(f"{a}:{co[a]}" for a in S if a in co)} }}  other={[a for a in co if a not in S]}')
