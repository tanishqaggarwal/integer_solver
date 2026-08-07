"""S11 step 91: release a GATED constant pin, then re-solve everything behind it.

§133-134 identified the circuit's residual as one two-condition primitive, and pin3.py
confirmed the closed form exactly:

    A = (w3-w1)^2*(w5+w1+w3+K) - (w4-w2)^2      matches x35389 to the digit
    B = (w6+w2)*(w3-w1) - (w1-w5)*(w4-w2)       matches x6671  to the digit

with ALL SEVEN quantities literal constants of the instance:

    w1 = x22152  w2 = x33462     pinned by a31670, a31672   -- gated by x24601
    w3 = x6418   w4 = x12553     pinned by a3576,  a3578    -- gated by x2081
    w5 = x22162  w6 = x30213     pinned by a1618,  a688     (through x24468, x18956)
    K  = x24453                  pinned by a41332           (bare, ungated)

On this branch the pair is therefore over-determined by its own constants and
does not close -- which is exactly why every repair in this lab has been conserved.
But four of the seven pins are GATED: `x24601*(w1 - C)` and `x2081*(w3 - C)`.  Set the
gate to zero and the value is free again.

x2081 and x4287 have been on this lab's FORBID list since Session 9 because flipping
them "cheats"; the price of x2081 = 0 was measured at 16 against a gain of 7.  That
was measured before the advice graph was solved and before A and B were
identified.  Re-price it here: release, let the handles absorb everything the release
leaves ≡ 0 mod p, then re-solve the freed values so A and B both vanish.

Usage: release.py [gate] [state.json]      gate in {2081, 24601, both}
"""
import os, sys
HERE = '/home/user/integer_solver/solve_lab/s10'
LAB = '/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
from intad import jacZ
import suppfree
P = ad.P
gate = sys.argv[1] if len(sys.argv) > 1 else '24601'
src = sys.argv[2] if len(sys.argv) > 2 else 'PIN_39013.json'
GATES = {'2081': [2081], '24601': [24601], 'both': [2081, 24601]}[gate]
v = L.load(src if os.path.isabs(src) else os.path.join(HERE, src))
ad.fwd(v, rounds=6)


def show(v, tag):
    av = L.all_atom_values(v)
    s = L.NEQ - len(L.failing_eqs(av))
    nz = [a for a in range(L.NA) if a not in L.atom_out and av[a]]
    print('%-32s score %-6d nonzero checks %-3d %s'
          % (tag, s, len(nz), nz[:14]), flush=True)
    return s, av, nz


show(v, 'start')
for g in GATES:
    v[g] = 0
ad.fwd(v, rounds=6)
s, av, nz = show(v, 'gates %s -> 0' % GATES)
print('   of those, ≡0 mod p (liftable): %s'
      % [a for a in nz if av[a] % P == 0], flush=True)

_, fl, S = suppfree.build(v, modp=None)
for it in range(40):
    av = L.all_atom_values(v)
    todo = [a for a in range(L.NA) if a not in L.atom_out and av[a] and av[a] % P == 0]
    if not todo:
        break
    cur = L.NEQ - len(L.failing_eqs(av))
    moved = False
    for c in todo:
        m = suppfree.atom_supp(c, v, S, modp=None)
        for i in range(len(fl)):
            if not ((m >> i) & 1):
                continue
            u = fl[i]
            g = jacZ(u, v, [c]).get(c, 0)
            if not g or g % P or av[c] % g:
                continue
            w = list(v)
            w[u] = w[u] - av[c] // g
            ad.fwd(w, rounds=6)
            a2 = L.all_atom_values(w)
            if a2[c] == 0 and L.NEQ - len(L.failing_eqs(a2)) >= cur:
                v, av, moved = w, a2, True
                break
        if moved:
            break
    if not moved:
        break
s, av, nz = show(v, 'after the integer lift')
T.save(v, os.path.join(HERE, 'REL_%s_%d.json' % (gate, s)))
print('saved REL_%s_%d.json' % (gate, s), flush=True)

# the freed values, and what A = B = 0 needs of them
FREEDBY = {2081: [(6418, 'w3'), (12553, 'w4')], 24601: [(22152, 'w1'), (33462, 'w2')]}
print('\ncoordinates released:')
for g in GATES:
    for t, nm in FREEDBY[g]:
        print('   x%-6d (%s) = %d' % (t, nm, v[t] % P))
print('\nA = x35389 ≡ %d\nB = x6671  ≡ %d' % (v[35389] % P, v[6671] % P), flush=True)


def probe(delta):
    w = list(v)
    for u, d in delta.items():
        w[u] = w[u] + d
    ad.fwd(w, rounds=6)
    return w, [w[35389] % P, w[6671] % P]


KN = [t for g in GATES for t, _ in FREEDBY[g]]
print('\nexact jacobian of (A, B) in the released knobs %s:' % KN, flush=True)
_, b = probe({})
cols = []
for u in KN:
    _, c = probe({u: 1})
    cols.append([(c[i] - b[i]) % P for i in range(2)])
M = [[cols[j][i] for j in range(len(KN))] for i in range(2)]
for r in M:
    print('   ', [str(x)[:26] for x in r])
if len(KN) >= 2:
    for j1 in range(len(KN)):
        for j2 in range(j1 + 1, len(KN)):
            det = (M[0][j1] * M[1][j2] - M[0][j2] * M[1][j1]) % P
            if not det:
                continue
            inv = pow(det, -1, P)
            r0, r1 = (-b[0]) % P, (-b[1]) % P
            d1 = (M[1][j2] * r0 - M[0][j2] * r1) % P * inv % P
            d2 = (M[0][j1] * r1 - M[1][j1] * r0) % P * inv % P
            w, ch = probe({KN[j1]: d1, KN[j2]: d2})
            aw = L.all_atom_values(w)
            sc = L.NEQ - len(L.failing_eqs(aw))
            print('   solve with (x%d, x%d): A=%d B=%d -> score %d'
                  % (KN[j1], KN[j2], ch[0], ch[1], sc), flush=True)
            if sc > s:
                T.save(w, os.path.join(HERE, 'RELS_%d.json' % sc))
                print('      *** saved RELS_%d.json' % sc, flush=True)
