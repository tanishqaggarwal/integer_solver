"""S10 step 41: the p-wire — the single structure the whole trapdoor rests on.

My census (S10 section 8) showed every one of the 1,249 solo handles has
granularity EXACTLY p.  The reason is that each handle enters as  wire * handle
with `wire` one of the ~220 variables pinned to p.  If the wire carried 1 instead
of p, every handle would have granularity 1 and BOTH congruences would dissolve.

So: what exactly forces the wire to p, and what would breaking it cost?
"""
import os, sys, collections, json
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L
import tools as T
import ad

P = ad.P
atom_out = L.atom_out
v = L.load(os.path.join(HERE, 'forward_state.json'))

WIRE = [u for u in range(L.NVARS) if v[u] == P]
print(f'variables equal to p: {len(WIRE)}')
print(f'  sample: {[f"x_{u}" for u in WIRE[:20]]}')
freew = [u for u in WIRE if u not in L.definer]
print(f'  of which FREE inputs: {len(freew)} -> {[f"x_{u}" for u in freew]}')

# which atoms define the wire members, and what is the root?
print('\n=== how each wire member is defined ===')
kinds = collections.Counter()
roots = []
for u in WIRE:
    d = L.definer.get(u)
    if d is None:
        kinds['FREE'] += 1; roots.append((u, None)); continue
    src = L.atom_src[d]
    others = [w for w in L.avars[d] if w != u]
    if not others:
        kinds['bare constant pin'] += 1; roots.append((u, d))
    elif all(w in WIRE for w in others):
        kinds['copy of another wire member'] += 1
    else:
        kinds['mixed'] += 1
for k, n in kinds.most_common():
    print(f'   {k:<32} {n}')
print('\nroot pins (no other variable in the defining atom):')
for u, d in roots:
    if d is None:
        print(f'   x_{u}: FREE INPUT (value p by choice!)')
    else:
        print(f'   x_{u}: a{d} = {L.atom_src[d][:150]}   neq={len(L.atom2eq.get(d,{}))}')

# how much of the instance touches the wire?
watoms = set()
for u in WIRE:
    watoms |= set(L.var_atoms[u])
weqs = set()
for a in watoms:
    weqs |= set(L.atom2eq.get(a, ()))
print(f'\natoms touching the wire: {len(watoms)}   equations: {len(weqs)}')

# what breaks if we set the whole wire to 1?
print('\n=== hypothetical: set every wire member to 1 ===')
w1 = list(v)
for u in WIRE:
    w1[u] = 1
av = L.all_atom_values(w1)
nz = [a for a in range(L.NA) if av[a]]
print(f'   nonzero atoms immediately: {len(nz)}')
ad.fwd(w1, rounds=3)
av = L.all_atom_values(w1)
nz = [a for a in range(L.NA) if av[a]]
nzc = [a for a in nz if a not in atom_out]
fail = L.failing_eqs(av)
print(f'   after forward eval: nonzero atoms={len(nz)} (checks {len(nzc)}) '
      f'failing={len(fail)} score={L.NEQ-len(fail)}')
print(f'   surviving checks: {nzc[:30]}')
T.save(w1, os.path.join(HERE, 'wire1.json'))
