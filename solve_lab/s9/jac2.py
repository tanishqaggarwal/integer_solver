"""Corrected Jacobian: square checks represented by their degree-2 ROOT."""
import pickle, time, sys
import harness as H
exec(open('repair.py').read().split('if __name__')[0])
NV=38748
roots = pickle.load(open('roots.pkl','rb'))
checks = [a for a in range(len(polys)) if a not in atom_out]
resid_poly = {a: (roots[a] if a in roots else polys[a]) for a in checks}
freeinp = [x for x in range(NV) if x not in definer]

def resid_all(v):
    return {a: evalpoly(P, v) for a, P in resid_poly.items()}

if __name__ == '__main__':
    v0 = H.load_assignment('S0.json')
    base = resid_all(v0)
    nzb = {a:x for a,x in base.items() if x}
    print('nonzero residuals at S0 (root model):', nzb)
    J = {}; t0=time.time()
    for i, f in enumerate(freeinp):
        v = list(v0); ripple(v, {f: v0[f]+1})
        col = {}
        for a, P in resid_poly.items():
            nv = evalpoly(P, v)
            if nv != base[a]: col[a] = nv - base[a]
        if col: J[f] = col
        if i % 1500 == 0: print(f'{i}/{len(freeinp)} {time.time()-t0:.0f}s', file=sys.stderr)
    pickle.dump({'J':J, 'base':base}, open('jac2.pkl','wb'))
    print(f'columns={len(J)} nnz={sum(len(c) for c in J.values())} time={time.time()-t0:.0f}s')
