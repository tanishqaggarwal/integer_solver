"""S10 step 77: ALL region knobs, not just solo handles.

eighth.py defined "adjustable" as carrying a solo free handle.  That is exactly
why it missed a22231, which is adjustable through x_28730 -- a variable already
in the placement.  Correct definition: a variable u is a REGION KNOB if moving it
changes no equation outside the twelve.

Scan every variable occurring in any atom of the twelve equations, perturb it
(with the residual atoms blocked so we control them by hand), and record its exact
effective equation footprint.  Every knob whose footprint lies inside the twelve
is a free parameter the Part I model never had.
"""
import os, sys, collections, json, time
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L
import tools as T

P = 2**256 - 2**32 - 977
v = L.load(os.path.join(LAB, 'best', 'new_instance_partial_39026.json'))
av = L.all_atom_values(v)
NZ = [22229, 22230, 35758, 35759, 35760, 35761, 35762]
E = sorted(L.eqs_of_atoms(NZ))
ES = set(E)
BLOCK = set(NZ) | {22231}
print(f'region: {len(E)} equations {E}')

# every atom in those equations, and every variable in those atoms
region_atoms = set()
for e in E:
    region_atoms |= set(L.eq_atoms[e][2])
cand = set()
for a in region_atoms:
    cand |= L.avars[a]
cand = sorted(cand)
print(f'atoms in the region: {len(region_atoms)}; variables in them: {len(cand)}')

base_av = av
BASE_NZ = set(a for a in range(L.NA) if av[a])
t0 = time.time()
knobs = []
for u in cand:
    w = list(v)
    try:
        ch, _ = L.ripple(w, {u: v[u] + 1}, block=BLOCK)
    except Exception:
        continue
    touched = L.touched_atoms(w, base_av, ch)
    if not touched:
        continue
    eqs = set()
    for a in touched:
        eqs |= set(L.atom2eq.get(a, ()))
    outside = eqs - ES
    moved_region_atoms = sorted(a for a in touched if a in region_atoms)
    knobs.append((len(outside), u, sorted(touched)[:8], moved_region_atoms[:8]))
knobs.sort()
print(f'scanned {len(cand)} variables in {time.time()-t0:.0f}s\n')
print(f'{"outside":>8} {"var":>9}  atoms moved (first 8)')
free_knobs = []
for out, u, touched, mv in knobs[:30]:
    tag = ''
    if out == 0:
        tag = '  <== FREE REGION KNOB'
        free_knobs.append(u)
    print(f'{out:>8} x_{u:<8} {touched}{tag}')
print(f'\nvariables whose ENTIRE footprint lies inside the twelve: {len(free_knobs)}')
print(f'  {[f"x_{u}" for u in free_knobs]}')

# which region atoms can those knobs move?
movable = set()
for u in free_knobs:
    w = list(v)
    ch, _ = L.ripple(w, {u: v[u] + 1}, block=BLOCK)
    for a in L.touched_atoms(w, base_av, ch):
        if a in region_atoms:
            movable.add(a)
print(f'\nregion atoms reachable by free knobs: {sorted(movable)}')
print(f'  (Part I used only {NZ})')
json.dump({'free_knobs': free_knobs, 'movable_atoms': sorted(movable),
           'region_eqs': E},
          open(os.path.join(HERE, 'regionknobs.json'), 'w'))
