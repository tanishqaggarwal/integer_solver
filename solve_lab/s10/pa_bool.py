"""Investigate the ultra-low-footprint boolean cluster {a33516,a33517,a27821,a27822}."""
import os, sys, collections, json
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L
C=[33516,33517,27821,27822,36244,36243,36245]
E=set()
for a in C[:4]: E|=set(L.atom2eq[a])
print('cluster',C[:4],'E=',sorted(E),'|E|=',len(E))
for e in sorted(E):
    m,sq,co=L.eq_atoms[e]
    print(f' eq{e} size={len(co)} atoms={sorted(co)}')
print()
allat=set()
for e in E: allat|=set(L.eq_atoms[e][2])
for a in sorted(allat):
    print(f'  a{a} fp={len(L.atom2eq[a])} in_E={len(set(L.atom2eq[a])&E)} out={len(set(L.atom2eq[a])-E)} gate={L.atom_out.get(a)} {L.atom_src[a][:80]}')
print()
v=L.load(os.path.join(LAB,'best','new_instance_partial_39026.json'))
for u in (29570,33095,24267,27026):
    print(f'x_{u}: val={v[u]} definer={L.definer.get(u)} nat={len(L.var_atoms[u])}')
    for a in L.var_atoms[u]:
        print(f'     a{a} fp={len(L.atom2eq[a])} gate={L.atom_out.get(a)} {L.atom_src[a][:110]}')
