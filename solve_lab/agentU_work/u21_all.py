"""U21: sweep every merge slot, one representative leaf pair per slot, both src choices."""
import sys, time, pickle, json, collections
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentU_work')
import u20_sweep as S
import umodel as U, uscore as SC

SLOTS = U.SLOTS
print('merge slots (both subtrees live): %d ; nodes in tree: %d ; ROOT=%d' % (len(SLOTS), len(U.order), U.ROOT))
out = {}
t0 = time.time()
for k, beta in enumerate(SLOTS):
    la, lb = U.tree[beta]
    A = sorted(U.LIVELEAF[la]); B = sorted(U.LIVELEAF[lb])
    a, b = A[0], B[0]
    rec = {'beta': beta, 'depth': U.depth[beta], 'nI': len(A), 'nJ': len(B), 'a': a, 'b': b}
    for tag, src in (('a', a), ('b', b)):
        try:
            PA, PB = S.coincide(beta, a, b, src)
            rec['coin_' + tag] = (PA == PB and PA is not None)
            n, vv, sd = S.price(beta, a, b, src)
            rec['n_' + tag] = n
        except Exception as e:
            rec['n_' + tag] = None; rec['err_' + tag] = repr(e)[:120]
    out[beta] = rec
    if k % 25 == 0:
        print('  %3d/%d beta=%d d=%d |I|=%d |J|=%d -> %s / %s   (%.0fs)'
              % (k, len(SLOTS), beta, rec['depth'], rec['nI'], rec['nJ'],
                 rec.get('n_a'), rec.get('n_b'), time.time() - t0))
        sys.stdout.flush()
pickle.dump(out, open('u_slot1.pkl', 'wb'))
vals = [min([x for x in (r.get('n_a'), r.get('n_b')) if x is not None], default=None) for r in out.values()]
vals = [v for v in vals if v is not None]
print('\n=== %d slots priced in %.0fs ===' % (len(vals), time.time() - t0))
print('min %d  median %d  max %d' % (min(vals), sorted(vals)[len(vals)//2], max(vals)))
print('distribution:', sorted(collections.Counter(vals).items())[:40])
best = sorted(out.values(), key=lambda r: min([x for x in (r.get('n_a'), r.get('n_b')) if x is not None], default=10**9))[:15]
print('\ncheapest 15 slots:')
for r in best:
    print('  beta=%-6d d=%-2d |I|=%-3d |J|=%-3d a=%d b=%d -> %s / %s  coin=%s/%s'
          % (r['beta'], r['depth'], r['nI'], r['nJ'], r['a'], r['b'],
             r.get('n_a'), r.get('n_b'), r.get('coin_a'), r.get('coin_b')))
