"""WR step 13: close a22229 on the uniform w=1 branch (would give 39020)."""
import os, sys, math, collections, random
HERE = '/home/user/integer_solver/solve_lab/s10'
LAB = '/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
import wr_frame as W
P = ad.P
F = W.F_WIRE
v = L.load(os.path.join(HERE, 'wr_rep1_1_39011.json'))
av = L.all_atom_values(v)
print('nonzero:', [a for a in range(L.NA) if av[a]])
print('a22229 =', av[22229])
print('a22229 src:', L.atom_src[22229])
D = v[7068] - v[2099]
print(f'x_7068 - x_2099 = {D}  mod 7376877 = {D % 7376877}')
print(f'x_642 = {v[642]}  7376877*x_642 = {7376877*v[642]}')
print(f'   need x_642 = (x_7068-x_2099)/7376877 -> divisible? {D % 7376877 == 0}')
for u in (642, 2099, 7068, 17325):
    print(f'  x_{u}: free={u in F.FREE} atoms={L.var_atoms[u]}')

print('\n--- route A: move x_642 (and x_17325 to keep a35762) ---')
if D % 7376877 == 0:
    w = list(v); w[642] = D // 7376877; w[17325] = w[642]
    F.fwd(w, rounds=8)
    F.report(w, '  route A')
else:
    print('  not divisible')

print('\n--- route B: adjust x_7068 so the difference becomes divisible ---')
# x_7068 is detached-free; changing it perturbs 8 atoms
for a in L.var_atoms[7068]:
    print(f'   a{a} neq={len(L.atom2eq[a])} {L.atom_src[a][:80]}')

print('\n--- route C: adjust x_2099 through its definer ---')
d = L.definer.get(2099)
print(f'   definer of x_2099 = a{d}: {L.atom_src[d] if d is not None else None}')
if d is not None:
    for u in sorted(L.avars[d]):
        print(f'      x_{u} free={u in F.FREE} atoms={len(L.var_atoms[u])}')

print('\n--- brute force: try every single-variable move that zeroes a22229 ---')
best = None
for u in sorted(L.avars[22229]) + [2099]:
    tgt = T.solve_lin(22229, u, v)
    if tgt is None:
        print(f'   x_{u}: not linearly solvable (non-integer or squared)')
        continue
    w = list(v); w[u] = tgt
    F.fwd(w, rounds=8)
    av2, nz2, fail2, sc2 = F.report(w, f'   set x_{u}')
    if best is None or sc2 > best[0]:
        best = (sc2, u, w)
if best:
    print(f'best single move: x_{best[1]} -> score {best[0]}')
    T.save(best[2], os.path.join(HERE, f'wr_a22229_{best[0]}.json'))
