"""S11 step 83: CONSTRUCT the residual away instead of searching for it.

§123 decompiled the seven residual atoms.  In frame 2 five of the seven have a
DETACHED output variable, so they are checks whose defining value we may simply
write down, and the other two reduce to one divisibility each:

    a35759   x29854  = 5113045*x7075*x9118          a35758  x29854 = x1329*p
    a35761   x31864  = -x7075*x8731                 a35760  x31864 = x10903*p
    a22230   x28730  = x9413*p        (x17499 = p)
    a35762   x642    = x17325*p       (x28599 = p)
    a22229   x7068   = x2099 + 7376877*x642

x1329, x9118, x10903, x8731, x9413, x17325 are free and x7068, x28730, x29854,
x31864, x642 are the detached ones, so ALL SEVEN can be written to zero exactly, the
moment p divides x7075*x9118 and x7075*x8731.  Forcing p | x9118 and p | x8731 does
that (the alternative, x7075 ≡ 0, is the already-priced `A = 0` route).

Nothing here is a search or a linearisation: every value is an exact integer written
from the identity it has to satisfy.  What it costs is whatever else those free
inputs feed, which only the checker can say.

Usage: build7.py [mode]      mode in {p9118, azero, both}
"""
import os, sys, itertools
HERE = '/home/user/integer_solver/solve_lab/s10'
LAB = '/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
from frame2 import definer, ORDER, FREE, CHECKS, fwd
P = ad.P
SEVEN = [22229, 22230, 35758, 35759, 35760, 35761, 35762]
mode = sys.argv[1] if len(sys.argv) > 1 else 'p9118'
base = L.load(os.path.join(LAB, 'best', 'new_instance_partial_39026.json'))
fwd(base)


def report(v, tag):
    av = L.all_atom_values(v)
    s = L.NEQ - len(L.failing_eqs(av))
    nz = [a for a in range(L.NA) if av[a]]
    seven = [a for a in SEVEN if av[a]]
    print('%-34s score %-6d nonzero atoms %-4d  residual still nonzero %s'
          % (tag, s, len(nz), seven), flush=True)
    return s, av


report(base, 'start (frame 2)')


def construct(v, force9118=True, force8731=True, extra=()):
    v = list(v)
    if force9118:
        v[9118] -= v[9118] % P
    if force8731:
        v[8731] -= v[8731] % P
    for u, val in extra:
        v[u] = val
    fwd(v)
    x7075 = v[7075]
    n1 = 5113045 * x7075 * v[9118]
    n2 = -x7075 * v[8731]
    if n1 % P or n2 % P:
        return None, 'p does not divide (%s, %s)' % (n1 % P == 0, n2 % P == 0)
    v[29854] = n1
    v[1329] = n1 // P
    v[31864] = n2
    v[10903] = n2 // P
    v[28730] = v[9413] * P
    v[642] = v[17325] * P
    fwd(v)
    v[7068] = v[2099] + 7376877 * v[642]
    fwd(v)
    return v, 'ok'


v, msg = construct(base)
print('construct: %s' % msg, flush=True)
if v is not None:
    s, av = report(v, 'after constructing all seven')
    T.save(v, os.path.join(HERE, 'B7_%d.json' % s))
    print('saved B7_%d.json' % s, flush=True)
    fail = sorted(L.failing_eqs(av))
    print('failing equations now: %s' % fail[:40], flush=True)
    nzc = [a for a in CHECKS if av[a]]
    nzg = [a for a in L.atom_out if av[a] and a not in set(SEVEN)]
    print('nonzero checks %s' % nzc[:20], flush=True)
    print('broken gates   %s' % nzg[:20], flush=True)

# variants: also drive the other free operands to convenient values
print('\n--- variants ---', flush=True)
best = 39026
for f9, f87 in itertools.product((True, False), repeat=2):
    for kx in ((), ((9413, 0),), ((17325, 0),), ((9413, 0), (17325, 0))):
        w, msg = construct(base, f9, f87, kx)
        if w is None:
            continue
        av = L.all_atom_values(w)
        s = L.NEQ - len(L.failing_eqs(av))
        tag = 'p9118=%s p8731=%s extra=%s' % (f9, f87, kx)
        print('  %-46s score %d  residual %s'
              % (tag, s, [a for a in SEVEN if av[a]]), flush=True)
        if s > best:
            best = s
            T.save(w, os.path.join(HERE, 'B7v_%d.json' % s))
            print('    *** saved B7v_%d.json' % s, flush=True)
