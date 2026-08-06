"""Global (block-free) ripple census over all variables, recording the PRIM-atom
footprint (square atoms represented by their degree-2 root) and a reverse index
atom -> variables that can move it.  Used as a fast candidate filter."""
import pickle, sys, time, collections
import lib as L
import model as MD

which = sys.argv[1] if len(sys.argv) > 1 else '24'
v0 = L.load(L.BEST24 if which == '24' else L.BEST22)
MD.BASEP = [MD.prim_val(a, v0) for a in range(L.NA)]
BASEP = MD.BASEP

DELTA = 1234567891011
foot = {}
rev = collections.defaultdict(set)
t0 = time.time()
for x in range(L.NVARS):
    v = list(v0)
    ch, st = L.ripple(v, {x: v0[x] + DELTA})
    cand = set()
    for u in ch:
        cand.update(MD.prim_var_atoms[u])
    tou = frozenset(a for a in cand if MD.prim_val(a, v) != BASEP[a])
    foot[x] = tou
    for a in tou:
        rev[a].add(x)
    if x % 8000 == 0:
        print(f'  {x}/{L.NVARS} {time.time()-t0:.0f}s', file=sys.stderr)
print(f'census2 done {time.time()-t0:.0f}s')
pickle.dump({'foot': foot, 'rev': dict(rev), 'base': BASEP},
            open(f'census2_{which}.pkl', 'wb'))
sz = collections.Counter(len(f) for f in foot.values())
print('prim-atom footprint size histogram (<=12):', {k: sz[k] for k in sorted(sz) if k <= 12})
