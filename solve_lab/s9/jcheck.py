"""Check how many Jacobian columns are invalid (a gate atom left unsatisfied by the ripple)."""
import pickle, time, sys
import harness as H
exec(open('repair.py').read().split('if __name__')[0])
NV=38748
gate_atoms = [a for a in range(len(polys)) if a in atom_out]
checks = [a for a in range(len(polys)) if a not in atom_out]
freeinp = [x for x in range(NV) if x not in definer]
v0 = H.load_assignment('S0.json')
bad = []; t0=time.time()
for i,f in enumerate(freeinp):
    v = list(v0); ch,_ = ripple(v, {f: v0[f]+1})
    # only gates whose vars changed can break
    broken = 0
    touched=set()
    for u in ch: touched.update(var_atoms[u])
    for a in touched:
        if a in atom_out and evalpoly(polys[a], v) != 0: broken += 1
    if broken: bad.append((f, broken))
    if i%1000==0: print(f'{i} t={time.time()-t0:.0f}s bad={len(bad)}', file=sys.stderr)
print(f'columns with broken gates: {len(bad)} / {len(freeinp)}')
print('examples:', bad[:20])
pickle.dump(bad, open('badcols.pkl','wb'))
