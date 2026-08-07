"""CL: apply the new machinery (K1/K2 gauge shifts, generic handle absorption,
Newton moves) to the DELIVERED 39,026 witness, where +1 would actually count."""
import os, sys, json, collections, itertools
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad
import cl_engine as E
P = E.P
definer, atom_out = L.definer, L.atom_out
FREE = set(ad.FREE)

W = os.path.join(LAB,'best','new_instance_partial_39026.json')
v0 = L.load(W)
av0, nz0, S0, bad0 = E.stats(v0)
print(f'delivered witness: score {S0}  nonzero atoms {nz0}')
print(f'failing eqs: {sorted(bad0)}')
touched = set()
for a in nz0: touched |= set(L.atom2eq.get(a,{}))
print(f'equations touched by the nonzero atoms: {len(touched)} -> {sorted(touched)}')
for a in nz0:
    print(f'  a{a:<6} eqs={sorted(L.atom2eq.get(a,{}))}  resid mod p = {av0[a]%P!=0}')
    print(f'      {L.atom_src[a][:150]}')

K1FREE = [8778, 14623, 16742, 31339, 33462]
K2FREE = [14853, 22152, 22649]
print('\n=== K-class values in this witness ===')
for u in K1FREE+K2FREE:
    print(f'  x_{u}: {v0[u].bit_length()} bits, mod p = {str(v0[u]%P)[:30]}..')

def trial(shift, tag):
    v = list(v0)
    for u,d in shift.items(): v[u] = v[u] + d
    ad.fwd(v, rounds=10)
    for _ in range(2):
        for a in [b for b in range(L.NA) if b not in atom_out and L.evalpoly(L.polys[b], v)]:
            E.absorb(v, a)
    av, nz, s, bad = E.stats(v)
    print(f'{tag}: score {s} ({s-S0:+d})  nz({len(nz)})={nz[:12]}')
    if s > S0:
        T.save(v, os.path.join(HERE,'cl_best.json'))
        print('   *** SAVED cl_best.json')
    return v, s

# 1. gauge shifts by a few natural deltas
for name, tgt in [('to x_24548', (v0[24548]-v0[14623])%P)]:
    trial({u: tgt for u in K1FREE}, f'K1 gauge {name}')
    trial({14623: tgt}, f'x_14623 only {name}')

# 2. generic absorption alone
trial({}, 'absorb-only')

# 3. Newton moves on every nonzero check
for a in nz0:
    if a in atom_out: continue
    vm = [x % P for x in v0]
    g = ad.grad(a, vm)
    rows = []
    for u in sorted(g):
        if u in E.FORBID: continue
        w = E.newton(v0, a, u, vm=vm, g=g, av=av0)
        if w is None: continue
        av, nz, s, bad = E.stats(w)
        rows.append((s,u,sorted(nz)))
        if s > S0:
            T.save(w, os.path.join(HERE,'cl_best.json'))
            print(f'   *** IMPROVED to {s} via x_{u} on a{a}; SAVED')
    rows.sort(key=lambda t:-t[0])
    print(f'\na{a} newton ({len(rows)} moves, support {len(g)}):')
    for s,u,nz in rows[:8]:
        print(f'   x_{u:<6} -> {s} ({s-S0:+d})  nz={nz[:9]}')
