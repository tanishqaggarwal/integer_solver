"""Certificate-guided minimum-equation-cost relaxation.
Repeatedly solve the mod-p system; when inconsistent, relax the certificate row whose
equation set adds least to the running union. Reports the achievable failing-equation bound."""
import pickle, collections, sys, os, time
HERE=os.path.dirname(os.path.abspath(__file__)); S9=os.path.dirname(HERE)
sys.path.insert(0,S9); os.chdir(S9)
P = 2**256-2**32-977
d = pickle.load(open('chain2/jac24.pkl','rb')); J = d['J']; base = d['base']
a2e = pickle.load(open('atom2eq.pkl','rb'))

def rows_of(J, base, drop):
    rows = collections.defaultdict(dict)
    for f, col in J.items():
        for a, dv in col.items():
            if a in drop: continue
            x = dv % P
            if x: rows[a][f] = x
    rhs = {a: (-base.get(a,0)) % P for a in rows}
    return rows, rhs

def eliminate(rows, rhs):
    """Return (None, certificate_provenance) if inconsistent, else (pivots, None)."""
    order = sorted(rows, key=lambda a: len(rows[a]))
    pivots = {}
    for a in order:
        row = dict(rows[a]); r = rhs[a]; pv = {a: 1}
        while True:
            hit = next((c for c in row if c in pivots), None)
            if hit is None: break
            prow, prhs, ppv = pivots[hit]
            fac = row[hit]*pow(prow[hit], P-2, P) % P
            for c, val in prow.items():
                nv = (row.get(c,0) - fac*val) % P
                if nv: row[c] = nv
                elif c in row: del row[c]
            r = (r - fac*prhs) % P
            for k, val in ppv.items(): pv[k] = (pv.get(k,0) - fac*val) % P
        if not row:
            if r: return None, {k: c for k, c in pv.items() if c}
            continue
        c = min(row); pivots[c] = (row, r, pv)
    return pivots, None

if __name__ == '__main__':
    drop = set()
    cost_union = set(a2e.get(22229, []))     # chain 1 already relaxed (atom 22229)
    print(f'start: chain-1 relaxed, equation cost = {len(cost_union)}')
    for it in range(30):
        rows, rhs = rows_of(J, base, drop)
        t0 = time.time()
        piv, cert = eliminate(rows, rhs)
        if cert is None:
            print(f'\nCONSISTENT after relaxing {sorted(drop)}  ({time.time()-t0:.0f}s)')
            print(f'  upper bound on failing equations = {len(cost_union)}')
            pickle.dump({'drop':sorted(drop),'cost':sorted(cost_union)}, open('chain2/relaxset.pkl','wb'))
            break
        # choose the certificate row that adds fewest new equations
        cands = []
        for a in cert:
            if a in drop: continue
            add = len(cost_union | set(a2e.get(a, []))) - len(cost_union)
            cands.append((add, len(a2e.get(a, [])), a))
        if not cands:
            print('no relaxable row in certificate'); break
        cands.sort()
        add, tot, a = cands[0]
        drop.add(a); cost_union |= set(a2e.get(a, []))
        print(f'it {it}: cert spans {len(cert)} rows; relax atom {a} '
              f'(+{add} eqs, atom has {tot}) -> running cost {len(cost_union)}   [{time.time()-t0:.0f}s]')
        if len(cost_union) >= 11:
            print('  >= 11 equations: no longer beats the 39,022 baseline; stopping.')
            break
