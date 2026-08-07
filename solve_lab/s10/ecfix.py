"""S11 step 88: the last two conditions are LINEAR in two unconstrained advice values.

§131 left three primitives at the advice fixed point (39,013):

    x15298*x11150 ≡ 0,   x15298*x25739 ≡ 0,   x15298*x37758 ≡ 0   (mod p)

and unfolding the operands shows x3023 and x2287 are not free but gate-defined:

    x3023 = 6926539*x6671        x2287 = 8272701*x35389

so with A = x35389 and B = x6671 the three conditions are HOMOGENEOUS LINEAR in
(A, B):

    x11150 =  8646263*A + 1073965*B      x25739 = 10159099*A + 6926539*B
    x37758 =  8272701*A + 5921311*B

a 3x2 system of rank 2, whose only solution is A ≡ B ≡ 0 (mod p).  And A and B are
the two point-addition identities:

    A = x35389 = x29322^2*x33469 - x3558^2         (x2-x1)^2*(x3+x1+x2) = (y2-y1)^2
    B = x6671  = x27713*x29322 - x1326*x3558       (y3+y1)*(x2-x1) = (y2-y1)*(x1-x3)

with x1 = x12186, y1 = x16742, x2 = x14853, y2 = x24908, x3 = x22162, y3 = x30213.
x22162 and x30213 are advice values whose only pins (a30976, a30978) are GATED by
x15574, which is zero -- so they are unconstrained, and A is linear in x22162 while
B is linear in both.  Two linear equations, two free unknowns.

The coefficients are recovered by exact evaluation rather than by hand: the maps are
linear, so one probe per unknown gives the exact column.

Usage: ecfix.py [state.json]
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
src = sys.argv[1] if len(sys.argv) > 1 else 'AG_39013.json'
v0 = L.load(src if os.path.isabs(src) else os.path.join(HERE, src))
ad.fwd(v0, rounds=6)
av = L.all_atom_values(v0)
print('%s: score %d; nonzero checks %s'
      % (src, L.NEQ - len(L.failing_eqs(av)),
         [a for a in range(L.NA) if a not in L.atom_out and av[a]]), flush=True)

KN = [22162, 30213]                      # the two unconstrained advice values
TGT = [35389, 6671]                      # the two quantities that must vanish mod p


def probe(delta):
    v = list(v0)
    for u, d in delta.items():
        v[u] = v[u] + d
    ad.fwd(v, rounds=6)
    return v, [v[t] % P for t in TGT]


_, b = probe({})
print('base: x35389 = %d,  x6671 = %d  (mod p)' % (b[0], b[1]), flush=True)
M = []
for u in KN:
    _, c = probe({u: 1})
    M.append([(c[i] - b[i]) % P for i in range(2)])
M = [[M[j][i] for j in range(len(KN))] for i in range(2)]     # rows = targets
print('exact 2x2 jacobian (rows x35389/x6671, cols x22162/x30213):')
for r in M:
    print('   ', r)
det = (M[0][0] * M[1][1] - M[0][1] * M[1][0]) % P
print('det = %d  (invertible: %s)' % (det, det != 0), flush=True)
if det == 0:
    print('singular -- the two knobs cannot both be solved'); sys.exit()
inv = pow(det, -1, P)
rhs = [(-b[0]) % P, (-b[1]) % P]
d0 = (M[1][1] * rhs[0] - M[0][1] * rhs[1]) % P * inv % P
d1 = (M[0][0] * rhs[1] - M[1][0] * rhs[0]) % P * inv % P
print('solution: x22162 += %d,  x30213 += %d  (mod p)' % (d0, d1), flush=True)

v, chk = probe({KN[0]: d0, KN[1]: d1})
print('after the jump: x35389 = %d, x6671 = %d' % (chk[0], chk[1]), flush=True)
av = L.all_atom_values(v)
for t in (11150, 25739, 37758):
    print('   x%-6d mod p = %d' % (t, v[t] % P))
nz = [a for a in range(L.NA) if a not in L.atom_out and av[a]]
print('score %d; nonzero checks %s; of which ≡0 mod p: %s'
      % (L.NEQ - len(L.failing_eqs(av)), nz, [a for a in nz if av[a] % P == 0]),
      flush=True)

# integer lift: every check that is now ≡0 mod p gets absorbed by a handle
_, freelist, SVS = suppfree.build(v, modp=None)
for rounds in range(6):
    av = L.all_atom_values(v)
    todo = [a for a in range(L.NA) if a not in L.atom_out and av[a] and av[a] % P == 0]
    if not todo:
        break
    moved = False
    for c in todo:
        m = suppfree.atom_supp(c, v, SVS, modp=None)
        for i in range(len(freelist)):
            if not ((m >> i) & 1):
                continue
            u = freelist[i]
            g = jacZ(u, v, [c]).get(c, 0)
            if not g or g % P or av[c] % g:
                continue
            w = list(v)
            w[u] = w[u] - av[c] // g
            ad.fwd(w, rounds=6)
            aw = L.all_atom_values(w)
            if aw[c] == 0:
                v, av, moved = w, aw, True
                print('   lift: a%d absorbed by x%d' % (c, u), flush=True)
                break
        if moved:
            break
    if not moved:
        break
av = L.all_atom_values(v)
s = L.NEQ - len(L.failing_eqs(av))
nz = [a for a in range(L.NA) if a not in L.atom_out and av[a]]
gb = [a for a in L.atom_out if av[a]]
print('\nFINAL score %d; nonzero checks %s; broken gates %d' % (s, nz, len(gb)))
T.save(v, os.path.join(HERE, 'EC_%d.json' % s))
print('saved EC_%d.json' % s)
