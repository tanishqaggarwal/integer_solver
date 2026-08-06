"""S12 step 21: cheap ripple-based repair of the RECORD frame (no forward-eval,
which destroys it).  For every nonzero atom, try every variable (and every free
input one level up through its definer) with solve_lin + L.ripple, and score.
Checkpoints to ac_r26fix.jsonl as it goes."""
import os, sys, json, time
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad
P = ad.P; definer = L.definer; FREE = set(ad.FREE); FORBID = {2081, 4287}
OUT = os.path.join(HERE, 'ac_r26fix.jsonl')
done = set()
if os.path.exists(OUT):
    for ln in open(OUT):
        try: done.add(tuple(json.loads(ln)['key']))
        except Exception: pass
v0 = L.load(os.path.join(LAB,'best','new_instance_partial_39026.json'))
av0 = L.all_atom_values(v0)
NZ = [a for a in range(L.NA) if av0[a]]
base = L.NEQ - len(L.failing_eqs(av0))
print(f'base {base}  nonzero {NZ}', flush=True)
f = open(OUT, 'a')
best = (base, None)
t0 = time.time(); n = 0
for a in NZ:
    cands = []
    for u in sorted(L.avars[a]):
        if u in FORBID: continue
        nv = T.solve_lin(a, u, v0)
        if nv is not None and nv != v0[u]: cands.append((u, nv, 'direct'))
        d = definer.get(u)
        if d is None: continue
        vv = list(v0)
        t2 = T.solve_lin(a, u, v0)
        if t2 is None: continue
        vv[u] = t2
        for z in sorted(L.avars[d]):
            if z == u or z in FORBID: continue
            nz2 = T.solve_lin(d, z, vv)
            if nz2 is not None and nz2 != v0[z]: cands.append((z, nz2, f'via{u}'))
    for u, nv, how in cands:
        key = [a, u, str(nv)[:40], how]
        if tuple(key) in done: continue
        w = list(v0)
        try: L.ripple(w, {u: nv})
        except Exception: continue
        av = L.all_atom_values(w)
        sc = L.NEQ - len(L.failing_eqs(av))
        n += 1
        f.write(json.dumps({'key': key, 'score': sc,
                            'nz': [b for b in range(L.NA) if av[b]][:20]}) + '\n'); f.flush()
        if sc > best[0]:
            best = (sc, (a, u, how))
            T.save(w, os.path.join(HERE, f'ac_r26_{sc}.json'))
            print(f'  IMPROVED: a{a} via x_{u} ({how}) -> {sc}', flush=True)
            if sc > 39026: T.save(w, os.path.join(HERE, 'ac_best.json'))
    print(f'  atom a{a}: {len(cands)} moves tried ({time.time()-t0:.0f}s) best {best[0]}', flush=True)
print(f'\n{n} moves scored; BEST {best[0]} {best[1]}')
