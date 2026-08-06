"""S11 step 21: can any DEAD path be activated cheaply?

113 free inputs reach the cluster structurally but have derivative zero, because
they only ever multiply a variable that is 0.  Setting one nonzero costs whatever
it switches on.  Scan them: which activate with little or no damage, and does the
cluster's gradient support GROW afterwards (i.e. do we gain a new knob)?
"""
import os, sys, random, json
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
P = ad.P
definer, atom_out = L.definer, L.atom_out
FREE = set(ad.FREE)
random.seed(11)
v0 = L.load(os.path.join(HERE, 'mod9118_0.json'))
vm0 = [x % P for x in v0]
av0 = L.all_atom_values(v0)
base_nz = set(a for a in range(L.NA) if av0[a])
BAD = [21617, 29539]
supp0 = set()
for a in BAD: supp0 |= set(ad.grad(a, vm0))
print(f'base score {L.NEQ-len(L.failing_eqs(av0))}; cluster gradient support {len(supp0)}')

def struct_cone(seeds):
    cone, st = set(), list(seeds)
    while st:
        t = st.pop()
        if t in cone: continue
        cone.add(t)
        a = definer.get(t)
        if a is None: continue
        for w in L.avars[a]:
            if w != t: st.append(w)
    return cone
seeds = set()
for a in BAD: seeds |= set(L.avars[a])
dead = sorted({u for u in struct_cone(seeds) if u in FREE and v0[u] == 0} - supp0)
print(f'dead free inputs to scan: {len(dead)}', flush=True)

res = []
for i, u in enumerate(dead):
    for val in (1, P):
        v = list(v0); v[u] = val
        ad.fwd(v, rounds=6)
        av = L.all_atom_values(v)
        nz = set(a for a in range(L.NA) if av[a])
        new = nz - base_nz
        s = L.NEQ - len(L.failing_eqs(av))
        grew = 0
        if not new or len(new) <= 2:
            vm = [x % P for x in v]
            supp = set()
            for a in BAD: supp |= set(ad.grad(a, vm))
            grew = len(supp - supp0)
        res.append((s, u, 'p' if val == P else '1', len(new), grew))
    if i % 25 == 0:
        print(f'  {i}/{len(dead)}', flush=True)
res.sort(reverse=True)
print('\ntop activations (score, input, value, new nonzero atoms, new knobs):')
for r in res[:16]: print(f'  {r}')
gain = [r for r in res if r[4] > 0]
print(f'\nactivations that GREW the cluster gradient support: {len(gain)}')
for r in gain[:12]: print(f'  {r}')
