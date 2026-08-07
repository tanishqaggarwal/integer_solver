"""S11 step 40: SECOND-ORDER activation.

Single activations provably cannot reach the cluster (s10/activate.py): a dead
free input u only multiplies a variable w that is 0, so d/du = 0 until w != 0.
Find, for each dead u, the blocking w, then a free input z that makes w nonzero,
and test the PAIR (u, z): does the cluster's gradient support grow?
"""
import os, sys, collections, time, random
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
P = ad.P
definer, atom_out = L.definer, L.atom_out
FREE = set(ad.FREE)
random.seed(23)
v0 = L.load(os.path.join(HERE, 'mod9118_0.json'))
vm0 = [x % P for x in v0]
av0 = L.all_atom_values(v0)
BAD = [21617, 29539]
supp0 = set()
for a in BAD: supp0 |= set(ad.grad(a, vm0))
base_nz = set(a for a in range(L.NA) if av0[a])
print(f'cluster gradient support {len(supp0)}', flush=True)

def cone(seeds):
    c, st = set(), list(seeds)
    while st:
        t = st.pop()
        if t in c: continue
        c.add(t)
        a = definer.get(t)
        if a is None: continue
        for w in L.avars[a]:
            if w != t: st.append(w)
    return c
seeds = set()
for a in BAD: seeds |= set(L.avars[a])
CC = cone(seeds)
dead = [u for u in CC if u in FREE and v0[u] == 0 and u not in supp0]

# blocking variables: w with v[w] == 0 that multiplies a dead u inside the cone
block = collections.Counter()
for u in dead:
    for a in L.var_atoms[u]:
        for m, c in L.polys[a].items():
            if u in m:
                for z in m:
                    if z != u and v0[z] == 0: block[z] += 1
print(f'dead free inputs {len(dead)}; blocking variables {len(block)}', flush=True)
top = [w for w, _ in block.most_common(30)]
# free inputs that can make a blocking variable nonzero
drivers = {}
for w in top:
    fr = [z for z in cone([w]) if z in FREE]
    if fr: drivers[w] = fr[:6]
print(f'blocking variables with free drivers: {len(drivers)}', flush=True)

t0 = time.time()
grew = []
tested = 0
for w, ds in list(drivers.items())[:14]:
    for z in ds[:3]:
        v = list(v0); v[z] = random.randrange(1, 1 << 64)
        ad.fwd(v, rounds=6)
        if v[w] == 0: continue
        vm = [x % P for x in v]
        supp = set()
        for a in BAD: supp |= set(ad.grad(a, vm))
        new = supp - supp0
        av = L.all_atom_values(v)
        nz = set(a for a in range(L.NA) if av[a]) - base_nz
        tested += 1
        if new:
            grew.append((len(new), z, w, len(nz)))
            print(f'  x_{z} activates x_{w}: support +{len(new)} '
                  f'(new atoms broken {len(nz)})', flush=True)
print(f'\ntested {tested} activations ({time.time()-t0:.0f}s); '
      f'{len(grew)} grew the cluster support')
for g in sorted(grew, reverse=True)[:10]: print(f'  {g}')
