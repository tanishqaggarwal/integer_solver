"""Full single-bit scan under both decoders + record structure."""
import sys, time, json, pickle
sys.path.insert(0, '/home/user/integer_solver/solve_lab/s9/sa')
import lib

lib.init_base()
out = {}
t0 = time.time()
for i, b in enumerate(lib.bfree):
    row = {}
    for al in (False, True):
        n, f, v, nz = lib.score([b], alignment=al)
        row['A' if al else 'N'] = (n, sorted(nz))
    out[b] = row
    if i % 200 == 0:
        print(f'{i}/{len(lib.bfree)} {time.time()-t0:.0f}s', file=sys.stderr)
pickle.dump(out, open('sa/scan1.pkl', 'wb'))

base = 11
for key, name in (('N', 'no-align'), ('A', 'align')):
    rows = sorted(out.items(), key=lambda t: t[1][key][0])
    print(f'\n=== {name}: best single bits ===')
    for b, r in rows[:20]:
        print(f'  x_{b}: {r[key][0]} fails, nz={len(r[key][1])} {r[key][1][:8]}')
    import collections
    c = collections.Counter(r[key][0] for r in out.values())
    print('  score histogram:', sorted(c.items())[:15])
