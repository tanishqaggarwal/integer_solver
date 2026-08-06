"""S10 step 78: boolean flips evaluated in the WITNESS frame, with region analysis.

Previous bit scans ran from the forward-eval state (37 failing) or used the
unreliable ripple.  The deliverable lives in a different frame (7 failing), and a
flip there changes WHICH pins are live -- e.g. x_24601 = 0 makes load pin 31670
vacuous, which would free x_22152 -> x_14853 -> x_7068 and could kill congruence 1
outright.

For each boolean free input: flip in the witness, ripple with the residual atoms
blocked, and report the resulting nonzero atoms, region size and failing count.
Then enumerate free region knobs for the most promising branches.
"""
import os, sys, collections, json, time
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L
import tools as T

BOOL = set()

P = 2**256 - 2**32 - 977
for _a, _poly in enumerate(L.polys):
    _ks = list(_poly.items())
    if len(_ks) == 2:
        _sq = [m for m, c in _ks if len(m) == 2 and m[0] == m[1]]
        _li = [m for m, c in _ks if len(m) == 1]
        if _sq and _li and _sq[0][0] == _li[0][0]:
            BOOL.add(_li[0][0])
base = L.load(os.path.join(LAB, 'best', 'new_instance_partial_39026.json'))
bav = L.all_atom_values(base)
NZ0 = [a for a in range(L.NA) if bav[a]]
BLOCK = set(NZ0) | {22231}
f0 = L.failing_eqs(bav)
print(f'witness: nonzero atoms {NZ0}, failing {len(f0)}, score {L.NEQ-len(f0)}')

FREEB = sorted(u for u in range(L.NVARS) if u not in L.definer and u in BOOL)
print(f'boolean free inputs: {len(FREEB)}')

t0 = time.time()
res = []
for i, b in enumerate(FREEB):
    v = list(base)
    try:
        L.ripple(v, {b: 1 - v[b]}, block=BLOCK)
    except Exception:
        continue
    av = L.all_atom_values(v)
    nz = [a for a in range(L.NA) if av[a]]
    fail = L.failing_eqs(av)
    region = L.eqs_of_atoms(nz)
    res.append((len(fail), len(nz), len(region), b))
    if len(fail) <= len(f0) + 4:
        print(f'  x_{b:<7} failing={len(fail):<4} nz_atoms={len(nz):<3} '
              f'region={len(region):<4} atoms={nz[:10]}', flush=True)
    if i % 200 == 0:
        print(f'  ... {i}/{len(FREEB)} ({time.time()-t0:.0f}s)', flush=True)

res.sort()
print(f'\n=== best 20 flips in the witness frame ({time.time()-t0:.0f}s) ===')
for f, na, nr, b in res[:20]:
    slack = nr - f          # how many region equations are already satisfied
    print(f'  x_{b:<7} failing={f:<5} score={L.NEQ-f:<7} nz_atoms={na:<3} '
          f'region={nr:<4} satisfied_in_region={slack}')
json.dump([{'bit': b, 'failing': f, 'natoms': na, 'region': nr}
           for f, na, nr, b in res[:80]],
          open(os.path.join(HERE, 'bitwitness.json'), 'w'))

# the three structural control bits, explicitly
print('\n=== the structural control bits ===')
for b in (24601, 2081, 4287):
    hit = [r for r in res if r[3] == b]
    if hit:
        f, na, nr, _ = hit[0]
        print(f'  x_{b}: failing={f} score={L.NEQ-f} nz_atoms={na} region={nr}')
    else:
        print(f'  x_{b}: not a boolean free input (definer '
              f'{L.definer.get(b)})')
