"""S11 step 82: exact univariate scan of the ELEVEN knobs at the deliverable.

eqker2 shows the residual at 39,026 is reachable by only ELEVEN non-handle free
inputs, and the 139 holding equations pin all eleven to first order (rank 11 of 11).
That is a linear statement; §122 says linear statements here bound nothing.  With
only eleven knobs the exact question is affordable:

for each knob u, interpolate the EXACT polynomial of every one of the 146 touched
equation-combinations along u, then
    G = gcd over the 139 that hold      -> its nonzero roots are jumps that break nothing
    H = gcd(G, the 7 that fail)         -> its nonzero roots WIN outright

and separately try the two jumps the decompilation of §123 asks for directly:
    x9118 -> x9118 - (x9118 mod p)      makes p | x7075*x9118, satisfying (†)
    x8731 -> x8731 - (x8731 mod p)      makes p | x7075*x8731, satisfying (‡)

Usage: f2uni.py [K]
"""
import os, sys, time
HERE = '/home/user/integer_solver/solve_lab/s10'
LAB = '/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
import suppfree
import fpoly as UP
from frame2 import definer, ORDER, FREE, CHECKS, fwd
P = ad.P
K = int(sys.argv[1]) if len(sys.argv) > 1 else 9

base = L.load(os.path.join(LAB, 'best', 'new_instance_partial_39026.json'))
fwd(base)
bm = [x % P for x in base]
bav = L.all_atom_values(base)
BASE = L.NEQ - len(L.failing_eqs(bav))
BAD = [a for a in CHECKS if bav[a]]
FAILEQ = sorted(L.failing_eqs(bav))
print(f'frame 2 @ {BASE}; failing equations {FAILEQ}', flush=True)

_, freelist, SVS = suppfree.build(base, definer=definer, ORDER=ORDER, FREE=FREE,
                                  modp=None)
U = set()
for c in BAD:
    m = suppfree.atom_supp(c, base, SVS, modp=None)
    U |= {freelist[i] for i in range(len(freelist)) if (m >> i) & 1}
U = sorted(U)
print(f'{len(U)} free inputs reach the residual: {U}', flush=True)
EQS = sorted(set().union(*[set(L.atom2eq[a]) for a in BAD]))
for u in U:                       # widen to every equation any knob can touch
    pass


def comb(av, e):
    s = 0
    for a, c in L.eq_atoms[e][2].items():
        if av[a]:
            s += c * av[a]
    return s % P


def scan(step, label):
    """step(v, d) mutates v by d units of the direction."""
    xs = list(range(K + 3))
    vals, EQ = [], None
    for d in xs:
        v = list(base)
        step(v, d)
        fwd(v)
        av = L.all_atom_values(v)
        if EQ is None:
            EQ = sorted(set(e for a in range(L.NA) if av[a] or bav[a]
                            for e in L.atom2eq[a]))
        vals.append([comb(av, e) for e in EQ])
    vary = [i for i in range(len(EQ)) if len(set(x[i] for x in vals)) > 1]
    G, fails, degbad = None, [], 0
    for i in vary:
        f = UP.interp(xs[:K + 1], [vals[d][i] for d in xs[:K + 1]])
        if not all(sum(co * pow(x, e2, P) for e2, co in enumerate(f)) % P == vals[x][i]
                   for x in xs[K + 1:]):
            degbad += 1
            continue
        (fails if vals[0][i] else []).append(f) if vals[0][i] else None
        if vals[0][i] == 0:
            G = f if G is None else UP.pgcd(G, f)
        else:
            fails.append(f)
    H = G
    for f in fails:
        H = f if H is None else UP.pgcd(H, f)
    sr = [r for r in UP.roots(G) if r % P] if G and len(G) > 1 else []
    hr = [r for r in UP.roots(H) if r % P] if H and len(H) > 1 else []
    best, bv = BASE, None
    for r in (sr + hr)[:10]:
        v = list(base)
        step(v, r)
        fwd(v)
        s = L.NEQ - len(L.failing_eqs(L.all_atom_values(v)))
        if s > best:
            best, bv = s, v
    print('  %-26s eqs %-4d vary %-4d degfail %-3d Gdeg %-4s safe-roots %-3d '
          'win-roots %-3d best %d'
          % (label, len(EQ), len(vary), degbad, (len(G) - 1) if G else None,
             len(sr), len(hr), best), flush=True)
    if bv is not None:
        T.save(bv, os.path.join(HERE, 'F2U_%d_%s.json' % (best, label)))
        print('    *** saved F2U_%d_%s.json' % (best, label), flush=True)
    return best


print('\n--- exact univariate scan of each knob ---', flush=True)
t0 = time.time()
for u in U:
    scan(lambda v, d, u=u: v.__setitem__(u, v[u] + d), 'x%d' % u)
print('(%.0fs)' % (time.time() - t0), flush=True)

print('\n--- the two jumps the decompilation asks for ---', flush=True)
for u in (9118, 8731, 1329, 10903, 9413, 17325):
    if u >= L.NVARS:
        continue
    d = (-base[u]) % P
    v = list(base)
    v[u] = v[u] + d
    fwd(v)
    av = L.all_atom_values(v)
    s = L.NEQ - len(L.failing_eqs(av))
    print('  x%-6d -> x%d = 0 mod p:  score %d  (nonzero atoms %d)'
          % (u, u, s, sum(1 for a in range(L.NA) if av[a])), flush=True)
    if s > BASE:
        T.save(v, os.path.join(HERE, 'F2Z_%d_x%d.json' % (s, u)))
        print('    *** saved', flush=True)
