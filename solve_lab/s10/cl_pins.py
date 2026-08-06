"""CL: the pin gadgets  x_2081 * (x_u - C)  and what x_2081 = 0 would cost."""
import os, sys, json, collections
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad
import cl_engine as E
P = E.P
definer, atom_out = L.definer, L.atom_out
FREE = set(ad.FREE)
v0 = L.load(os.path.join(HERE,'mod9118_0.json'))
av0, nz0, S0, bad0 = E.stats(v0)

print('=== atoms that are multiples of x_2081 / x_4287 ===')
pins = collections.defaultdict(list)
for a in range(L.NA):
    if a in atom_out: continue
    Pp = L.polys[a]
    for b in (2081, 4287):
        if all(b in m for m in Pp):
            pins[b].append(a)
for b, lst in pins.items():
    print(f'  x_{b}: {len(lst)} CHECK atoms are multiples of it, total eqs '
          f'{len(set().union(*[set(L.atom2eq.get(a,{})) for a in lst]))}')
    for a in lst[:40]:
        print(f'     a{a:<6} eqs={len(L.atom2eq.get(a,{})):<3} {L.atom_src[a][:120]}')

print('\n=== a31672 and a3576 ===')
for a in (31672, 3576, 25676, 42245, 33796):
    print(f'a{a}: eqs={len(L.atom2eq.get(a,{}))} {L.atom_src[a][:220]}')

print('\n=== how many equations involve x_2081 at all ===')
print(f'  var_eqs[2081] = {len(L.var_eqs[2081])}, var_eqs[4287] = {len(L.var_eqs[4287])}')

print('\n=== branch tests ===')
for b1 in (0,1):
    for b2 in (0,1):
        v = list(v0); v[2081]=b1; v[4287]=b2
        ad.fwd(v, rounds=10)
        av = L.all_atom_values(v)
        nz = [a for a in range(L.NA) if av[a]]
        bad = set(L.failing_eqs(av))
        s = L.NEQ-len(bad)
        print(f'  x_2081={b1} x_4287={b2}: score {s} ({s-S0:+d})  nz({len(nz)})={nz[:16]}')
