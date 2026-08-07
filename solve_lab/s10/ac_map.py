"""S12 step 3: the FULL activation map.

For every dead free input u in the cluster's structural cone, enumerate the
monomial-level BLOCKING SETS (the zero variables that kill u's derivative), then
for every blocking variable w the set of free inputs z that structurally reach w.
Saves ac_map.json for the pair sweep.
"""
import os, sys, json, time, collections
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad
P = ad.P
definer, atom_out = L.definer, L.atom_out
FREE = set(ad.FREE)
v0 = L.load(os.path.join(HERE,'mod9118_0.json'))
vm0 = [x % P for x in v0]
av0 = L.all_atom_values(v0)
BAD = [21617, 29539]
supp0 = set()
for a in BAD: supp0 |= set(ad.grad(a, vm0))
print(f'cluster gradient support {len(supp0)} free inputs', flush=True)

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
print(f'cluster structural cone: {len(CC)} vars', flush=True)
dead = sorted(u for u in CC if u in FREE and v0[u] == 0 and u not in supp0)
print(f'dead free inputs in the cone: {len(dead)}', flush=True)

# ---- monomial-level blocking sets --------------------------------------------
# atoms that matter: those defining a variable of the cone, or BAD themselves
CCatoms = set(BAD)
for t in CC:
    a = definer.get(t)
    if a is not None: CCatoms.add(a)
blocks = {}                       # u -> list of frozenset(blocking vars)
allW = collections.Counter()
for u in dead:
    bs = set()
    for a in L.var_atoms[u]:
        if a not in CCatoms: continue
        for m in L.polys[a]:
            if u not in m: continue
            z = frozenset(w for w in m if w != u and v0[w] == 0)
            if z: bs.add(z)
            else: bs.add(frozenset())   # already live path (shouldn't happen)
    if bs:
        blocks[u] = sorted(bs, key=len)
        for s in bs:
            for w in s: allW[w] += 1
print(f'dead inputs with blocking sets: {len(blocks)}; distinct blockers {len(allW)}', flush=True)
sz = collections.Counter(len(s) for u in blocks for s in blocks[u])
print(f'blocking-set sizes: {dict(sorted(sz.items()))}', flush=True)

# ---- which free inputs can drive each blocker --------------------------------
t0 = time.time()
W = sorted(allW)
drivers = {}
conecache = {}
for i, w in enumerate(W):
    c = cone([w])
    fr = sorted(z for z in c if z in FREE and v0[z] == 0)
    drivers[w] = fr
    if i % 100 == 0: print(f'  cone {i}/{len(W)} ({time.time()-t0:.0f}s)', flush=True)
pool = sorted(set().union(*[set(d) for d in drivers.values()])) if drivers else []
nd = collections.Counter(len(d) for d in drivers.values())
print(f'blockers with >=1 free driver: {sum(1 for d in drivers.values() if d)}/{len(W)}')
print(f'DRIVER POOL (distinct free inputs able to light some blocker): {len(pool)}')
print(f'driver-count histogram (blockers by #drivers): '
      f'{sorted(nd.items())[:6]} ... max {max(nd) if nd else 0}')
json.dump({'dead': dead, 'supp0': sorted(supp0), 'W': W,
           'blocks': {str(u): [sorted(s) for s in blocks[u]] for u in blocks},
           'drivers': {str(w): drivers[w] for w in W},
           'pool': pool},
          open(os.path.join(HERE,'ac_map.json'),'w'))
print(f'saved ac_map.json ({time.time()-t0:.0f}s)')
