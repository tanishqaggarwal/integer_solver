"""S10 step 7: the MUX branch, examined in EQUATION space.

On b1=b2=1 the MUX routes x_2099 <- x_9118 and x_19964 <- x_8731, both FREE, which
makes the two surviving congruences satisfiable.  Prior sessions rejected this branch
because it activates loads 19088/22233/22235 and those atoms cannot all vanish.
But atoms need not vanish -- EQUATIONS must.  Check the overlap.
"""
import os, sys, collections
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
import lib as L

P = 2**256 - 2**32 - 977
v = L.load(os.path.join(LAB, 'best', 'new_instance_partial_39026.json'))
av = L.all_atom_values(v)

print('=== how the MUX selector is built ===')
for u in (4287, 2081, 21279, 7075, 25297, 37158, 2099, 19964, 4432, 12553, 6418,
          31861, 14865, 9106, 2239, 31731, 4306, 27177):
    d = L.definer.get(u)
    src = L.atom_src[d][:120] if d is not None else 'FREE INPUT'
    print(f'  x_{u:<6} = {str(v[u])[:34]:<36} def_atom={str(d):<7} {src}')

LOADS = [19088, 22233, 22235, 19087, 22232, 22234]
RESID = [22229, 22230, 35758, 35759, 35760, 35761, 35762]
print('\n=== equations of the load atoms vs equations of the residual atoms ===')
lo = {a: set(L.atom2eq.get(a, {})) for a in LOADS}
re_ = {a: set(L.atom2eq.get(a, {})) for a in RESID}
allload = set().union(*lo.values())
allres = set().union(*re_.values())
for a in LOADS:
    print(f'  load a{a:<6} eqs={sorted(lo[a])}')
print(f'  union(load eqs)     = {sorted(allload)}  ({len(allload)})')
print(f'  union(residual eqs) = {sorted(allres)}  ({len(allres)})')
print(f'  load eqs NOT touched by residual atoms = {sorted(allload-allres)}')
print(f'  -> these are the equations the MUX branch would break with no residual atom')
print(f'     available to compensate.')

print('\n=== atoms of each equation in the union (how much room to compensate) ===')
for i in sorted(allload | allres):
    m, sq, co = L.eq_atoms[i]
    nzt = [a for a in co if av[a]]
    print(f'  eq {i:<6} n_atoms={len(co):<3} sq={int(sq)} nonzero_now={nzt}')
