"""S11 step 118: solve the divisibility lock, then absorb with HANDLES ONLY.

divlock shows the last lock is solvable: the step of `(10159099*a + 6926539*b) mod
6672769` per unit of k is 1963712 for x22162 and 3063958 for x30213, and
gcd(those, 6672769) = 1, so the congruence has solutions.

But every solved candidate still scored 39,015 or less, because the generic absorber
is greedy and habsorb had already found the trap: **x22162 absorbs a1618 with
coefficient 1, and x30213 absorbs a688 with coefficient 8863713**.  Those are the
values themselves, so "absorbing" a1618 that way just puts w5 back on its pin and
destroys A = 0.  The absorber has to be restricted to genuine HANDLES -- free inputs
whose exact integer coefficient on the check is a multiple of p -- which is what makes
the move invisible mod p and therefore harmless to every congruence.

Order: shift (k1, k2) to solve the lock, then absorb a19297 (x30317), a19299 (x5146)
and a30984 (x2936), then anything else a handle can take.

Usage: finish.py [state.json] [K1LO] [K1HI]
"""
import os, sys, math, time
HERE = '/home/user/integer_solver/solve_lab/s10'
LAB = '/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
from intad import jacZ
import suppfree
P = ad.P
Q = 6672769
src = sys.argv[1] if len(sys.argv) > 1 else 'PF_best_39015.json'
K1LO = int(sys.argv[2]) if len(sys.argv) > 2 else -30
K1HI = int(sys.argv[3]) if len(sys.argv) > 3 else 30
base = L.load(src if os.path.isabs(src) else os.path.join(HERE, src))
ad.fwd(base, rounds=6)


def report(v, tag):
    av = L.all_atom_values(v)
    f = L.failing_eqs(av)
    s = L.NEQ - len(f)
    nz = [a for a in range(L.NA) if a not in L.atom_out and av[a]]
    print('%-42s score %-6d failing %-3d checks %s' % (tag, s, len(f), nz),
          flush=True)
    return s, av, nz


BASE, av0, NZ = report(base, 'start (%s)' % os.path.basename(src))
_, freelist, SVS = suppfree.build(base, modp=None)
FORBID = {22162, 30213, 14853, 16742, 12186, 24908}


def handle_absorb(v, rounds=25):
    """Absorb only through genuine handles: coefficient a multiple of p."""
    for _ in range(rounds):
        av = L.all_atom_values(v)
        cur = L.NEQ - len(L.failing_eqs(av))
        moved = False
        for c in [x for x in range(L.NA) if x not in L.atom_out and av[x]]:
            m = suppfree.atom_supp(c, v, SVS, modp=None)
            for i in range(len(freelist)):
                if not ((m >> i) & 1):
                    continue
                u = freelist[i]
                if u in FORBID:
                    continue
                col = jacZ(u, v, [c]).get(c, 0)
                if not col or col % P or av[c] % col:
                    continue
                w = list(v)
                w[u] = w[u] - av[c] // col
                ad.fwd(w, rounds=6)
                aw = L.all_atom_values(w)
                if aw[c] == 0 and L.NEQ - len(L.failing_eqs(aw)) >= cur:
                    v, moved = w, True
                    break
            if moved:
                break
        if not moved:
            break
    return v


print('\n--- control: handle-only absorption at the base state ---', flush=True)
v = handle_absorb(list(base))
best, bestv = report(v, 'handles only, no shift')[0], list(v)
if best > BASE:
    T.save(v, os.path.join(HERE, 'FIN_%d.json' % best))
    print('   *** NEW BEST %d' % best, flush=True)

A, B = base[35389], base[6671]
if A % P or B % P:
    print('A or B not ≡ 0 mod p at the base'); sys.exit()
a, b = A // P, B // P
tgt = (10159099 * a + 6926539 * b) % Q
EFF = {}
for u in (22162, 30213):
    w = list(base)
    w[u] += P
    ad.fwd(w, rounds=6)
    EFF[u] = ((w[35389] - A) // P, (w[6671] - B) // P)
c1 = (10159099 * EFF[22162][0] + 6926539 * EFF[22162][1]) % Q
c2 = (10159099 * EFF[30213][0] + 6926539 * EFF[30213][1]) % Q
print('\nlock target %d; steps %d and %d (mod %d)' % (tgt, c1, c2, Q), flush=True)

print('\n--- shift to solve the lock, then absorb with handles only ---', flush=True)
t0 = time.time()
tried = 0
for n1 in range(K1LO, K1HI + 1):
    if time.time() - t0 > 1800:
        break
    r = (tgt + c1 * n1) % Q
    g = math.gcd(c2, Q)
    if r % g:
        continue
    n2 = (-(r // g)) * pow(c2 // g, -1, Q // g) % (Q // g)
    for n2v in (n2, n2 - Q // g):
        v = list(base)
        v[22162] += n1 * P
        v[30213] += n2v * P
        ad.fwd(v, rounds=6)
        if v[35389] % P or v[6671] % P:
            continue
        lock = (10159099 * (v[35389] // P) + 6926539 * (v[6671] // P)) % Q
        v = handle_absorb(v)
        s, aw, nz = report(v, '   k1=%d k2=%s  lock %d'
                           % (n1, str(n2v)[:14], lock))
        tried += 1
        if s > best:
            best, bestv = s, list(v)
            T.save(v, os.path.join(HERE, 'FIN_%d.json' % s))
            print('      *** NEW BEST %d -- saved FIN_%d.json' % (s, s), flush=True)
print('\n%d shifts tried; best %d (was %d)' % (tried, best, BASE))
