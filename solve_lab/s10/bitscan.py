"""S10 step 40: systematic boolean scan in the CORRECT forward-eval frame.

Every prior bit scan used the ripple (which I showed is unreliable).  Here each
flip is followed by a genuine topological forward evaluation, so the measured
state is a real circuit state.  A bit flip changes WHICH pins are live, i.e. the
structure of the system -- the only freedom left after the linearisation was
shown to be rigid.
"""
import os, sys, json, time, collections
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L
import tools as T
import ad
from newton import BOOL

P = ad.P
atom_out = L.atom_out
base = L.load(os.path.join(HERE, 'forward_state.json'))
BFREE = sorted(u for u in ad.FREE if u in BOOL)
print(f'boolean free inputs: {len(BFREE)}', flush=True)

av0 = L.all_atom_values(base)
f0 = L.failing_eqs(av0)
print(f'base: failing={len(f0)} score={L.NEQ-len(f0)} '
      f'checks={[a for a in range(L.NA) if av0[a] and a not in atom_out]}', flush=True)

t0 = time.time()
res = []
for i, b in enumerate(BFREE):
    v = list(base)
    v[b] = 1 - v[b]
    ad.fwd(v, rounds=2)
    av = L.all_atom_values(v)
    nzc = [a for a in range(L.NA) if av[a] and a not in atom_out]
    nzg = [a for a in range(L.NA) if av[a] and a in atom_out]
    fail = L.failing_eqs(av)
    res.append((len(fail), len(nzc), b, nzc[:12], len(nzg)))
    if len(fail) <= len(f0):
        print(f'  x_{b}: failing={len(fail)} score={L.NEQ-len(fail)} '
              f'checks={nzc[:12]} gates_broken={len(nzg)}', flush=True)
    if i % 100 == 0:
        print(f'  ... {i}/{len(BFREE)}  ({time.time()-t0:.0f}s)', flush=True)

res.sort()
print(f'\n=== best 25 single flips ({time.time()-t0:.0f}s) ===', flush=True)
for f, nc, b, chk, ng in res[:25]:
    print(f'  x_{b:<7} failing={f:<5} score={L.NEQ-f:<7} nchecks={nc:<3} '
          f'gates_broken={ng:<3} {chk}')
json.dump([{'bit': b, 'failing': f, 'nchecks': nc, 'checks': chk}
           for f, nc, b, chk, ng in res[:200]],
          open(os.path.join(HERE, 'bitscan.json'), 'w'))
print('saved bitscan.json')
