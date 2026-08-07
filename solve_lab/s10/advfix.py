"""S11 step 86: set the two wrong advice values and let the handles absorb.

§85: the free content of the instance is THIRTEEN 296-bit advice values, each of
which must satisfy one congruence x_i ≡ y_i (mod p) enforced by a gadget
c*(x_i - y_i) - p*h.  At B7_39004 eleven of the thirteen hold; exactly two fail:

    a7930   9367949 *(x24548 - x25442) - x7927     ->  x24548 ≡ x25442 (mod p)
    a29539  12846437*(x14853 - x1308 ) - x29967    ->  x14853 ≡ x1308  (mod p)

and the two remaining nonzero checks, a41512 and a40826, contain exactly those two
residuals (a41512 carries -27*a7930 - 16*x7927, and x7927 is p times a handle), so
they should follow for free.

Each advice value is k*p + r.  Six of the thirteen share k = 839192594282 and three
share k = 1094785891323, so k is structured, not random -- which matters because the
congruence fixes only r, leaving k for us to choose.  Try the congruence-correct r
against every k that occurs in the instance.

Usage: advfix.py [state.json]
"""
import os, sys, itertools
HERE = '/home/user/integer_solver/solve_lab/s10'
LAB = '/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
from intad import jacZ
import suppfree
P = ad.P
src = sys.argv[1] if len(sys.argv) > 1 else 'B7_39004.json'
v0 = L.load(src if os.path.isabs(src) else os.path.join(HERE, src))
ad.fwd(v0, rounds=6)
FREE = [t for t in range(L.NVARS) if t not in L.definer]
ADV = [t for t in FREE if 280 <= v0[t].bit_length() <= 330]
KS = sorted({v0[t] // P for t in ADV})
print('%s: score %d;  ks in use %s'
      % (src, L.NEQ - len(L.failing_eqs(L.all_atom_values(v0))), KS), flush=True)
_, freelist, SVS = suppfree.build(v0, modp=None)


def handles(c, v):
    m = suppfree.atom_supp(c, v, SVS, modp=None)
    out = []
    for i in range(len(freelist)):
        if (m >> i) & 1:
            u = freelist[i]
            g = jacZ(u, v, [c]).get(c, 0)
            if g and g % P == 0:
                out.append((u, g))
    return out


# (check, advice input, target variable, coefficient)
GAD = [(7930, 24548, 25442, 9367949), (29539, 14853, 1308, 12846437)]


def apply(v, picks):
    v = list(v)
    for (c, x, y, co), k in picks:
        v[x] = k * P + (v[y] % P)
    ad.fwd(v, rounds=6)
    av = L.all_atom_values(v)
    for c, x, y, co in [g for g, _ in picks]:
        if av[c] % P:
            continue
        for h, g in handles(c, v):
            if av[c] % g == 0:
                v[h] = v[h] - av[c] // g
                ad.fwd(v, rounds=6)
                av = L.all_atom_values(v)
                break
    return v, av


best = L.NEQ - len(L.failing_eqs(L.all_atom_values(v0)))
print('\n--- one at a time ---', flush=True)
for g in GAD:
    for k in KS + [v0[g[1]] // P]:
        v, av = apply(v0, [(g, k)])
        s = L.NEQ - len(L.failing_eqs(av))
        nz = [a for a in range(L.NA) if a not in L.atom_out and av[a]]
        gb = sum(1 for a in L.atom_out if av[a])
        print('  a%-6d x%-6d k=%-14d score %-6d checks %s gates %d'
              % (g[0], g[1], k, s, nz[:8], gb), flush=True)
        if s > best:
            best = s
            T.save(v, os.path.join(HERE, 'ADV_%d.json' % s))
            print('    *** saved ADV_%d.json' % s, flush=True)

print('\n--- both together ---', flush=True)
for k1 in KS:
    for k2 in KS:
        v, av = apply(v0, [(GAD[0], k1), (GAD[1], k2)])
        s = L.NEQ - len(L.failing_eqs(av))
        nz = [a for a in range(L.NA) if a not in L.atom_out and av[a]]
        print('  k1=%-14d k2=%-14d score %-6d checks %s'
              % (k1, k2, s, nz[:8]), flush=True)
        if s > best:
            best = s
            T.save(v, os.path.join(HERE, 'ADV2_%d.json' % s))
            print('    *** saved ADV2_%d.json' % s, flush=True)
print('\nbest %d' % best)
