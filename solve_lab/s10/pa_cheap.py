import os, sys, collections, json
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L
# primitives = atoms with footprint>1 ; bundles = footprint 1
bund=[a for a in range(L.NA) if len(L.atom2eq[a])==1]
prim=[a for a in range(L.NA) if len(L.atom2eq[a])>1]
print('bundles',len(bund),'primitives',len(prim))
# verify bundle equations have size 1
bad=[b for b in bund if len(L.eq_atoms[next(iter(L.atom2eq[b]))][2])!=1]
print('bundles whose equation has >1 atom:',len(bad))
# cheap primitives
cheap=sorted(prim,key=lambda a:len(L.atom2eq[a]))[:30]
for a in cheap:
    print(f'a{a} neq={len(L.atom2eq[a])} gate={L.atom_out.get(a)} nvars={len(L.avars[a])} eqs={sorted(L.atom2eq[a])}')
    print(f'    {L.atom_src[a][:150]}')
# free variable count
print('\nfree vars:', L.NVARS-len(L.definer), 'definers', len(L.definer))
v = L.load(os.path.join(LAB,'best','new_instance_partial_39026.json'))
print('nonzero vars in witness:', sum(1 for x in v if x!=0))
print('x_26064 =', v[26064])
P=2**256-2**32-977
print('P=',P)
