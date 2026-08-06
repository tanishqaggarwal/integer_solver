"""Newton solve mod p at a point where the core gradient is non-degenerate."""
import pickle, time, sys, collections
import harness as H
import modsolve
exec(open('repair.py').read().split('if __name__')[0])
P = 2**256-2**32-977
NV = 38748
roots = pickle.load(open('roots.pkl','rb'))
checks = [a for a in range(len(polys)) if a not in atom_out]
resid_poly = {a: (roots[a] if a in roots else polys[a]) for a in checks}
freeinp = [x for x in range(NV) if x not in definer]

def resid(v):
    return {a: evalpoly(Pp, v) for a, Pp in resid_poly.items()}

def jac_at(v0, base, cols=None):
    J = {}
    for f in (cols if cols is not None else freeinp):
        v = list(v0); ripple(v, {f: v0[f]+1})
        col = {}
        for a, Pp in resid_poly.items():
            nv = evalpoly(Pp, v)
            if nv != base[a]: col[a] = nv - base[a]
        if col: J[f] = col
    return J

def build_S1():
    v = H.load_assignment('S0.json')
    # fix both mirror congruences by moving the free leaves
    ripple(v, {14853: v[1308], 24548: v[25442]})
    return v

if __name__ == '__main__':
    v = build_S1()
    H.save_assignment(v, 'S1.json')
    base = resid(v)
    print('S1 nonzero residuals:', {a:x for a,x in base.items() if x})
    t0=time.time()
    J = jac_at(v, base)
    print(f'jacobian: {len(J)} cols, nnz={sum(len(c) for c in J.values())}, {time.time()-t0:.0f}s')
    pickle.dump({'J':J,'base':base}, open('jacS1.pkl','wb'))
    sol, msg = modsolve.solve(J, base)
    print('solve:', msg)
    if sol is not None:
        nz={c:x for c,x in sol.items() if x}
        print(f'{len(nz)} nonzero deltas mod p')
        pickle.dump(nz, open('modsolS1.pkl','wb'))
