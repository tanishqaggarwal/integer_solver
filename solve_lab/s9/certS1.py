"""Find the certificate of inconsistency: which rows combine to 0 = nonzero."""
import pickle, collections, time
P = 2**256-2**32-977
d = pickle.load(open('jacS1.pkl','rb')); J = d['J']; base = d['base']
rows = collections.defaultdict(dict)
for f, col in J.items():
    for a, dv in col.items():
        x = dv % P
        if x: rows[a][f] = x
rhs = {a: (-base.get(a, 0)) % P for a in rows}
order = sorted(rows, key=lambda a: len(rows[a]))
pivots = {}
prov = {}     # col -> {orig_row: coeff}
for a in order:
    row = dict(rows[a]); r = rhs[a]; pv = {a: 1}
    while True:
        hit = None
        for c in row:
            if c in pivots: hit = c; break
        if hit is None: break
        prow, prhs, ppv = pivots[hit]
        factor = row[hit]*pow(prow[hit], P-2, P) % P
        for c, val in prow.items():
            nv = (row.get(c,0) - factor*val) % P
            if nv: row[c] = nv
            elif c in row: del row[c]
        r = (r - factor*prhs) % P
        for k, val in ppv.items(): pv[k] = (pv.get(k,0) - factor*val) % P
    if not row:
        if r:
            print(f'CONTRADICTION reached while processing atom {a}, residual {r}')
            nzprov = {k: c for k, c in pv.items() if c}
            print(f'  involves {len(nzprov)} original rows:')
            for k, c in sorted(nzprov.items()):
                print(f'    atom {k}  coeff={c}   base_residual={base.get(k,0)}  (nonzero base? {base.get(k,0)!=0})')
            break
        continue
    c = min(row)
    pivots[c] = (row, r, pv)
