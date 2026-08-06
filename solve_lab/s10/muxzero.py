"""S10 step 22: on the MUX branch, zero ALL SEVEN residual atoms (possible there,
since x_9118 drives x_2099 and x_8731 drives x_19964) and see what is left."""
import os, sys, json
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L
import tools as T

P = 2**256 - 2**32 - 977
NZ = [22229, 22230, 35758, 35759, 35760, 35761, 35762]
BLOCK = set(NZ) | {22231}
v = L.load(os.path.join(HERE, 'mux_on.json'))
print('x_7075 =', v[7075], ' x_21279 =', v[21279])

# close the two load pins x_4287 lit, using their p-handles
for a, u, h in ((3568, 31861, 6504), (3570, 14865, 26658)):
    c, rest = T.lin_parts(a, u, v)
    nv = -rest // c
    L.ripple(v, {u: nv}, block=BLOCK)
    print(f'  a{a}: set x_{u} -> value  (atom now {L.evalpoly(L.polys[a], v) if hasattr(L,"evalpoly") else "?"})')

av = L.all_atom_values(v)
print('after pinning MUX inputs: nz =', [a for a in range(L.NA) if av[a]])

# now zero the seven residual atoms
seeds = {29854: 0, 1329: 0, 31864: 0, 10903: 0, 642: 0, 17325: 0}
seeds[9118] = v[7068] - v[37158]        # -> x_2099 = x_7068  => a22229 = 0
seeds[8731] = v[4432] - v[20492]        # -> x_19964 = x_4432 => x_28730 = 0
seeds[28730] = 0
seeds[9413] = 0
L.ripple(v, seeds, block=BLOCK)
av = L.all_atom_values(v)
nz = [a for a in range(L.NA) if av[a]]
fail = L.failing_eqs(av)
print(f'\nresidual atoms now: {[(a, av[a]) for a in NZ]}')
print(f'nonzero atoms: {nz}')
print(f'failing: {len(fail)} -> score {L.NEQ-len(fail)}')
print(f'failing eqs: {fail[:40]}')
T.save(v, os.path.join(HERE, 'muxzero.json'))
