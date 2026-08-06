"""S11 step 113: the coarse generators cost exactly ONE advice congruence each.

repaircost pins down what moving the two coarse generators actually breaks:

    x28730 += 1  ->  16 equations, whose only nonzero atoms are a7930 and a41512
                     -- i.e. the advice congruence  x24548 = x25442
    x7068  += 1  ->  13 equations, whose only nonzero atoms are a29539 and a40826
                     -- i.e. the advice congruence  x14853 = x1308

and x24548 and x14853 are FREE advice values, so each congruence is repaired by a
residue jump plus a handle -- the two-phase move of §125, which is exact.  If that
works, a0 and a1 become free generators, the alpha-lattice is complete, and every one
of the twelve rows becomes reachable.

    a22229 = x7068 - x2099 - 7376877*x642        a0 fine via x7068
    a22230 = x28730 - x9413*x17499               a1 fine via x28730
    a7930  = 9367949 *(x24548 - x25442) - x7927  x7927  = p*handle
    a29539 = 12846437*(x14853 - x1308 ) - x29967 x29967 = p*handle

Usage: coarse.py [d7068] [d28730]
"""
import os, sys, time
HERE = '/home/user/integer_solver/solve_lab/s10'
LAB = '/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
from frame2 import definer, ORDER, FREE, fwd
from intad import jacZ
import suppfree
P = ad.P
d7068 = int(sys.argv[1]) if len(sys.argv) > 1 else 1
d28730 = int(sys.argv[2]) if len(sys.argv) > 2 else 0
SEVEN = [22229, 22230, 35758, 35759, 35760, 35761, 35762]
base = L.load(os.path.join(LAB, 'best', 'new_instance_partial_39026.json'))
fwd(base)
E12 = set(e for a in SEVEN for e in L.atom2eq[a])


def report(v, tag):
    av = L.all_atom_values(v)
    f = set(L.failing_eqs(av))
    s = L.NEQ - len(f)
    print('%-38s score %-6d  inside-12 %-3d  outside %-3d  %s'
          % (tag, s, len(f & E12), len(f - E12), sorted(f - E12)[:10]), flush=True)
    return s, av


report(base, 'witness')
v = list(base)
v[7068] += d7068
v[28730] += d28730
fwd(v)
report(v, 'after the coarse move')

# repair each broken advice congruence: residue jump on the free advice value,
# then the handle absorbs the quotient over Z
GAD = [(7930, 24548, 25442, 9367949), (29539, 14853, 1308, 12846437)]
_, freelist, SVS = suppfree.build(v, definer=definer, ORDER=ORDER, FREE=FREE,
                                  modp=None)


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


for c, x, y, co in GAD:
    av = L.all_atom_values(v)
    if av[c] == 0:
        continue
    d = (v[x] - v[y]) % P
    v[x] -= d
    fwd(v)
    av = L.all_atom_values(v)
    if av[c] % P == 0 and av[c]:
        for h, g in handles(c, v):
            if av[c] % g == 0:
                v[h] = v[h] - av[c] // g
                fwd(v)
                av = L.all_atom_values(v)
                break
    report(v, 'after repairing a%d via x%d' % (c, x))

# the repair moves an advice value, and the advice DAG has to follow: every
# downstream congruence x_i = y_i must be re-imposed, in topological order.
TWO = [(8778, 16144), (24548, 25442), (14623, 27522), (16742, 19083),
       (22649, 29524), (31339, 6858), (14853, 1308)]
for sweep_i in range(4):
    changed = 0
    for t, y in TWO:
        d = (v[t] - v[y]) % P
        if d:
            v[t] -= d
            changed += 1
            fwd(v)
    s, av = report(v, 'advice sweep %d (changed %d)' % (sweep_i, changed))
    if not changed:
        break
# integer lift for anything now ≡ 0 (mod p)
for _ in range(20):
    av = L.all_atom_values(v)
    todo = [a for a in range(L.NA) if a not in L.atom_out and av[a]
            and av[a] % P == 0]
    cur = L.NEQ - len(L.failing_eqs(av))
    moved = False
    for c in todo:
        for h, g in handles(c, v):
            if av[c] % g == 0:
                w = list(v)
                w[h] = w[h] - av[c] // g
                fwd(w)
                a2 = L.all_atom_values(w)
                if a2[c] == 0 and L.NEQ - len(L.failing_eqs(a2)) >= cur:
                    v, moved = w, True
                    break
        if moved:
            break
    if not moved:
        break
s, av = report(v, 'FINAL')
print('\nthe seven residual atoms now:')
for a in SEVEN:
    print('   a%-6d %s' % (a, 'ZERO' if av[a] == 0 else str(av[a])[:34] + '...'))
T.save(v, os.path.join(HERE, 'CO_%d_%d_%d.json' % (d7068, d28730, s)))
print('saved CO_%d_%d_%d.json' % (d7068, d28730, s))
