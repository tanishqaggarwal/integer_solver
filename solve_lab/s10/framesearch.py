"""S11 step 60: search the FRAME space itself.

Every measurement so far fixed a frame (canonical, frame2, frame3, F_WIRE) and
searched inside it.  The break census sampled single-variable perturbations.  This
samples the frames: pick a random set of variables to DETACH (so their defining
atoms become checks free to be nonzero), forward-evaluate from the deliverable with
a block-preserving pass, run the enriched engine, record the score.

Chunked and resumable via chunk.sweep -- usage: framesearch.py START END [SEED]
"""
import os, sys, random, time
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad
from chunk import sweep
P = ad.P
SEVEN = [22229, 22230, 35758, 35759, 35760, 35761, 35762]
BASE_DET = {7068: 22229, 28730: 22230, 29854: 35758, 31864: 35761, 642: 35762}
base = L.load(os.path.join(LAB, 'best', 'new_instance_partial_39026.json'))
av0 = L.all_atom_values(base)
# candidate variables to additionally detach: outputs of gate atoms near the residual
cone = set()
st = list(SEVEN)
seen = set()
while st:
    a = st.pop()
    if a in seen: continue
    seen.add(a)
    for w in L.avars[a]:
        cone.add(w)
        for b in L.var_atoms[w]:
            if b not in seen and len(seen) < 400: st.append(b)
CAND = sorted(w for w in cone if w in L.definer)
print(f'candidate variables to detach: {len(CAND)}', flush=True)

def evaluate(spec):
    extra = spec
    DET = dict(BASE_DET)
    for w in extra:
        a = L.definer.get(w)
        if a is not None: DET[w] = a
    definer = {t: a for t, a in L.definer.items() if t not in DET}
    ORDER = [t for t in ad.ORDER if t not in DET]
    v = list(base)
    for _ in range(4):
        for u in ORDER:
            nv = T.solve_lin(definer[u], u, v)
            if nv is not None: v[u] = nv
    av = L.all_atom_values(v)
    f = len(L.failing_eqs(av))
    nz = [a for a in range(L.NA) if av[a]]
    rec = {'score': L.NEQ - f, 'nz': len(nz), 'det': list(extra)}
    if L.NEQ - f > 39026:
        T.save(v, os.path.join(HERE, f'FS_{L.NEQ-f}.json'))
        rec['SAVED'] = True
    return rec

start = int(sys.argv[1]); end = int(sys.argv[2])
seed = int(sys.argv[3]) if len(sys.argv) > 3 else 0
random.seed(seed * 100003 + start)
specs = []
for i in range(end - start):
    k = random.randint(1, 5)
    specs.append(tuple(sorted(random.sample(CAND, min(k, len(CAND))))))
sweep(f'framesearch_s{seed}', specs, evaluate, 0, len(specs),
      keyfn=lambda s: ','.join(map(str, s)), budget=500)
best = max((r['score'] for r in __import__('chunk').load(f'framesearch_s{seed}')),
           default=0)
print(f'best score this seed so far: {best}')
