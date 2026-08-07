import os, sys, collections, json
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad

print('NA', L.NA, 'NEQ', L.NEQ, 'NVARS', L.NVARS)
gates = set(L.atom_out)
checks = [a for a in range(L.NA) if a not in gates]
print('gate atoms', len(gates), 'check atoms', len(checks))
fp = collections.Counter()
for a in range(L.NA):
    fp[len(L.atom2eq.get(a,{}))]+=1
print('footprint histogram (all atoms):', sorted(fp.items())[:15])
fpc = collections.Counter(len(L.atom2eq.get(a,{})) for a in checks)
print('footprint histogram (checks):', sorted(fpc.items())[:15])
# atoms not in any equation
print('atoms with 0 equations:', sum(1 for a in range(L.NA) if not L.atom2eq.get(a)))

v = L.load(os.path.join(LAB,'best','new_instance_partial_39026.json'))
av = L.all_atom_values(v)
nz = [a for a in range(L.NA) if av[a]]
print('best-state nonzero atoms:', nz)
print('failing eqs:', L.failing_eqs(av))
for a in nz:
    print(f'  a{a} gate={a in gates} out={L.atom_out.get(a)} neq={len(L.atom2eq.get(a,{}))} eqs={sorted(L.atom2eq.get(a,{}))}')
