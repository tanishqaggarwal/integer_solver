"""Sensitivity (exact affine Jacobian) at S0 for all free inputs."""
import pickle, time, sys
import harness as H
exec(open('repair.py').read().split('if __name__')[0])
P = 2**256-2**32-977
NV = 38748
checks = [a for a in range(len(polys)) if a not in atom_out]
freeinp = [x for x in range(NV) if x not in definer]

if __name__ == '__main__':
    v0 = H.load_assignment('S0.json')
    base = {a: evalpoly(polys[a], v0) for a in checks}
    J = {}
    t0 = time.time()
    for i, f in enumerate(freeinp):
        v = list(v0)
        ripple(v, {f: v0[f] + 1})
        col = {}
        for a in checks:
            nv = evalpoly(polys[a], v)
            if nv != base[a]: col[a] = nv - base[a]
        if col: J[f] = col
        if i % 1000 == 0: print(f'{i}/{len(freeinp)} {time.time()-t0:.0f}s', file=sys.stderr)
    pickle.dump({'J':J, 'base':base}, open('jac.pkl','wb'))
    print(f'free inputs with nonzero column: {len(J)}; nnz={sum(len(c) for c in J.values())}')
    print(f'time {time.time()-t0:.0f}s')
