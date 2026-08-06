"""S11 step 9: are the VETO (absorbable) rows linear mod p under large moves?

If yes, the reduced closure is not a linearisation -- it is exact, and this
stratum is closed by linear algebra.  If no, the veto is untrustworthy and a
large move can escape it.  Also: show the actual collateral of a pair solve.
"""
import os, sys, random, collections, json
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
from fwdad import jac_column
P = ad.P
definer, atom_out = L.definer, L.atom_out
FREE = set(ad.FREE); FORBID = {2081, 4287}
random.seed(7)
v0 = L.load(os.path.join(HERE, 'mod9118_0.json'))
vm0 = [x % P for x in v0]
av0 = L.all_atom_values(v0)
CHECKS = sorted(a for a in range(L.NA) if a not in atom_out)
ABS = set(json.load(open(os.path.join(HERE, 'absorbable.json')))['absorbable'])
BAD = [21617, 29539]
U = sorted((set(ad.grad(BAD[0], vm0)) | set(ad.grad(BAD[1], vm0))) - FORBID,
           key=lambda u: len(L.var_atoms[u]))

print('=== linearity of veto rows under LARGE moves ===')
ok = bad = 0
badrows = collections.Counter()
for u in U[:16]:
    col = jac_column(u, v0, vm0, CHECKS)
    d = random.randrange(1, P)
    w = list(v0); w[u] = w[u] + d
    ad.fwd(w, rounds=8)
    aw = L.all_atom_values(w)
    for c in col:
        if c not in ABS: continue
        pred = (av0[c] + col[c] * d) % P
        if pred == aw[c] % P: ok += 1
        else: bad += 1; badrows[c] += 1
print(f'absorbable rows: linear {ok}, NONLINEAR {bad} '
      f'({100*bad/max(1,ok+bad):.1f}%)')
print(f'  rows that broke linearity: {len(badrows)} {list(badrows)[:10]}')

print('\n=== collateral of one exact pair solve ===')
g = {a: ad.grad(a, vm0) for a in BAD}
r = [(-av0[a]) % P for a in BAD]
done = 0
for i in range(len(U)):
    for j in range(i + 1, len(U)):
        u1, u2 = U[i], U[j]
        a11, a21 = g[BAD[0]].get(u1, 0) % P, g[BAD[1]].get(u1, 0) % P
        a12, a22 = g[BAD[0]].get(u2, 0) % P, g[BAD[1]].get(u2, 0) % P
        det = (a11 * a22 - a12 * a21) % P
        if det == 0: continue
        di = pow(det, -1, P)
        d1 = (a22 * r[0] - a12 * r[1]) % P * di % P
        d2 = (a11 * r[1] - a21 * r[0]) % P * di % P
        w = list(v0); w[u1] += d1; w[u2] += d2
        ad.fwd(w, rounds=8)
        aw = L.all_atom_values(w)
        nz = [a for a in range(L.NA) if aw[a]]
        nzabs = [a for a in nz if a in ABS]
        print(f'  x_{u1}+x_{u2}: a21617 {aw[21617]%P==0} a29539 {aw[29539]%P==0}  '
              f'nonzero {len(nz)} (absorbable {len(nzabs)})  '
              f'failing {len(L.failing_eqs(aw))}')
        done += 1
        if done >= 6: break
    if done >= 6: break
