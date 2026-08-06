"""S11 step 90: x3 and y3 are pinned to LITERAL CONSTANTS -- the circuit's output.

Unfolding a688 and a1618 with the selectors at their current values:

    x32237 = x21023*x22820,  x21023 = p, x22820 free   -> a handle, ≡ 0 (mod p)
    x34243 = x14393*x16153,  x16153 = p, x14393 free   -> a handle, ≡ 0 (mod p)
    x25538 = x16742*x34606 + x5647*x24908,  x34606 = x5647 = 0   -> 0
    x13913 = x12186*x34606 + x5647*x14853                        -> 0

so, with x15298 = 1,

    x18956 ≡ y3 (mod p)        x24468 ≡ x3 (mod p)

and the two pins read off directly:

    a688   8863713*x18956 ≡ C1   ->   y3 ≡ C1 * 8863713^{-1}  (mod p)
    a1618  x24468        ≡ C2   ->   x3 ≡ C2                 (mod p)

Both C1 and C2 are literals in EQUATIONS.txt.  So x3 and y3 are not free after all --
they are the circuit's OUTPUT POINT, written into the instance.  §134 solved A = B = 0
by moving them, which is why the pins broke.  Set them to the pinned values instead
and the remaining question becomes the honest one: does the point addition close for
the input coordinates the rest of the circuit forces?

Usage: pin3.py [state.json]
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
C1 = 1114942656963403660822546820446916783439088877768247923308647546252105232931473698035897478439338
C2 = 91416258160755509149180373473728639746431157665678710450404458852172057265575180278101002
src = sys.argv[1] if len(sys.argv) > 1 else 'EC_39014.json'
v = L.load(src if os.path.isabs(src) else os.path.join(HERE, src))
ad.fwd(v, rounds=6)
av = L.all_atom_values(v)
print('%s: score %d; nonzero checks %s'
      % (src, L.NEQ - len(L.failing_eqs(av)),
         [a for a in range(L.NA) if a not in L.atom_out and av[a]]), flush=True)
y3req = C1 % P * pow(8863713, -1, P) % P
x3req = C2 % P
print('pinned  x3 = x22162 ≡ %d' % x3req)
print('pinned  y3 = x30213 ≡ %d' % y3req, flush=True)
print('current x22162 ≡ %d  (match %s)' % (v[22162] % P, v[22162] % P == x3req))
print('current x30213 ≡ %d  (match %s)' % (v[30213] % P, v[30213] % P == y3req),
      flush=True)

k1, k2 = v[22162] // P, v[30213] // P
v[22162] = k1 * P + x3req
v[30213] = k2 * P + y3req
ad.fwd(v, rounds=6)
av = L.all_atom_values(v)
print('\nafter pinning: x18956 ≡ %d, x24468 ≡ %d' % (v[18956] % P, v[24468] % P))
print('   a688 = %d,  a1618 = %d' % (av[688] % P, av[1618] % P))
print('   A = x35389 ≡ %d' % (v[35389] % P))
print('   B = x6671  ≡ %d' % (v[6671] % P), flush=True)

_, fl, S = suppfree.build(v, modp=None)
for _ in range(10):
    av = L.all_atom_values(v)
    todo = [a for a in range(L.NA) if a not in L.atom_out and av[a] and av[a] % P == 0]
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
            if a2[c] == 0 and L.NEQ - len(L.failing_eqs(a2)) >= \
                              L.NEQ - len(L.failing_eqs(av)):
                v, av, moved = w, a2, True
                print('   lift: a%d absorbed by x%d' % (c, u), flush=True)
                break
        if moved:
            break
    if not moved:
        break
av = L.all_atom_values(v)
s = L.NEQ - len(L.failing_eqs(av))
nz = [a for a in range(L.NA) if a not in L.atom_out and av[a]]
print('\nFINAL score %d; nonzero checks %s; broken gates %d'
      % (s, nz, sum(1 for a in L.atom_out if av[a])))
T.save(v, os.path.join(HERE, 'PIN_%d.json' % s))
print('saved PIN_%d.json' % s)
