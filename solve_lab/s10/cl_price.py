"""CL step 3: price moving the RHS free inputs x_14623 / x_14853 (+ handle absorption)."""
import os, sys, json, collections
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad
P = 2**256 - 2**32 - 977
definer, atom_out = L.definer, L.atom_out
FREE = set(ad.FREE)

v0 = L.load(os.path.join(HERE,'mod9118_0.json'))
av0 = L.all_atom_values(v0)
def score(v):
    return L.NEQ - len(L.failing_eqs(L.all_atom_values(v)))
S0 = score(v0)
print(f'base frame score {S0}, nonzero atoms {[a for a in range(L.NA) if av0[a]]}')

print('\n--- handle variables ---')
for h in (986, 5040, 11360, 30163):
    print(f'x_{h}: {"FREE" if h in FREE else "gate a"+str(definer[h])}  val={str(v0[h])[:40]} bits={v0[h].bit_length()} ==p? {v0[h]==P} consumers={len(L.var_atoms[h])} -> {sorted(L.var_atoms[h])[:10]}')

# ---------------- price: fix a21617 by moving x_14623 --------------------
def try_rhs(target, gate, rhs, handle, wire):
    print(f'\n===== fix a{target} by moving free input x_{rhs}, absorb into x_{handle} =====')
    r = av0[target] % P
    vm = [x % P for x in v0]
    g = ad.grad(target, vm)
    d = g.get(rhs)
    print(f'  residue mod p = {r}')
    print(f'  d(a{target})/d(x_{rhs}) = {d}')
    if d is None: return None
    delta = (-r * pow(d,-1,P)) % P
    v = list(v0); v[rhs] = v[rhs] + delta
    ad.fwd(v, rounds=8)
    av = L.all_atom_values(v)
    print(f'  after move+fwd: a{target} mod p = {av[target]%P}   (should be 0)')
    # absorb: solve check for the gate output var, then set the handle
    w = atom_out[gate][1]
    tgt = T.solve_lin(target, w, v)
    print(f'  need x_{w} = {str(tgt)[:40]}...  (currently {v[w]})')
    if tgt is not None:
        vv = list(v); vv[w] = tgt
        nv = T.solve_lin(gate, handle, vv)
        if nv is None:
            print(f'  handle x_{handle} cannot realise it')
        else:
            v[handle] = nv
            ad.fwd(v, rounds=8)
    av = L.all_atom_values(v)
    nz = [a for a in range(L.NA) if av[a]]
    s = L.NEQ - len(L.failing_eqs(av))
    bad0 = set(L.failing_eqs(av0)); bad1 = set(L.failing_eqs(av))
    print(f'  SCORE {s} ({s-S0:+d})   nonzero atoms {nz}')
    print(f'  eqs fixed: {len(bad0-bad1)}   eqs broken: {len(bad1-bad0)}')
    print(f'  newly nonzero atoms: {sorted(set(nz)-set(a for a in range(L.NA) if av0[a]))[:40]}')
    return v, s

r1 = try_rhs(21617, 21616, 14623, 5040, 986)
r2 = try_rhs(29539, 29538, 14853, 30163, 11360)

# ---------------- both at once --------------------
if r1 and r2:
    print('\n===== BOTH RHS moved =====')
    vm = [x % P for x in v0]
    v = list(v0)
    for target, gate, rhs, handle in [(21617,21616,14623,5040),(29539,29538,14853,30163)]:
        g = ad.grad(target, vm); r = av0[target] % P
        d = g[rhs]; v[rhs] = v[rhs] + (-r*pow(d,-1,P))%P
    ad.fwd(v, rounds=8)
    av = L.all_atom_values(v)
    print(f'  residues mod p: a21617 {av[21617]%P}  a29539 {av[29539]%P}')
    for target, gate, handle in [(21617,21616,5040),(29539,29538,30163)]:
        w = atom_out[gate][1]
        tgt = T.solve_lin(target, w, v)
        if tgt is None: continue
        vv = list(v); vv[w]=tgt
        nv = T.solve_lin(gate, handle, vv)
        if nv is not None: v[handle]=nv
    ad.fwd(v, rounds=8)
    av = L.all_atom_values(v)
    nz=[a for a in range(L.NA) if av[a]]
    s = L.NEQ-len(L.failing_eqs(av))
    print(f'  SCORE {s} ({s-S0:+d})  nonzero {nz[:40]} (n={len(nz)})')
    T.save(v, os.path.join(HERE,'cl_bothrhs.json'))
