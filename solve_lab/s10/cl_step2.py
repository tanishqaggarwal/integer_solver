"""CL: apply the zero-cost x_6418 repair of a29539, then look at what appeared (a3576)."""
import os, sys, json, time
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad
import cl_engine as E
P = E.P

v0 = L.load(os.path.join(HERE,'mod9118_0.json'))
av0, nz0, S0, bad0 = E.stats(v0)
print(f'base: score {S0} nz {nz0}')

v1 = E.newton(v0, 29539, 6418)
av1, nz1, S1, bad1 = E.stats(v1)
print(f'after x_6418 fix of a29539: score {S1} nz {nz1}  broke {len(bad1-bad0)} fixed {len(bad0-bad1)}')
T.save(v1, os.path.join(HERE,'cl_s1.json'))

for a in nz1:
    print(f'\na{a}: check={E.is_check(a)} eqs={len(L.atom2eq.get(a,{}))} {sorted(L.atom2eq.get(a,{}))}')
    print(f'   src: {L.atom_src[a][:180]}')
    print(f'   handles: {[(f"x_{w}",f"a{d}",f"wire x_{wi}",f"h x_{h}") for w,d,wi,h in E.handles_of(a, v1)]}')
    vm = [x % P for x in v1]
    g = ad.grad(a, vm)
    print(f'   grad support {len(g)}')

# now price every Newton move on a3576 and a21617 from state v1
for tgt in nz1:
    if not E.is_check(tgt): continue
    vm = [x % P for x in v1]
    g = ad.grad(tgt, vm)
    rows = []
    for u in sorted(g):
        if u in E.FORBID: continue
        w = E.newton(v1, tgt, u, vm=vm, g=g, av=av1)
        if w is None: continue
        av, nz, s, bad = E.stats(w)
        rows.append((s, u, sorted(nz), len(bad-bad1), len(bad1-bad)))
    rows.sort(key=lambda t: -t[0])
    print(f'\n=== from cl_s1, repairing a{tgt} ({len(rows)} moves) ===')
    for s,u,nz,br,fx in rows[:12]:
        print(f'  x_{u:<6} -> score {s:>6} ({s-S1:+d})  broke {br} fixed {fx}  nz={nz[:10]}')
    json.dump([[r[0],r[1],r[2],r[3],r[4]] for r in rows], open(os.path.join(HERE,f'cl_s1_{tgt}.json'),'w'))
