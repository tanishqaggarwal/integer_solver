"""S11 step 92: close the addition THROUGH the advice solve.

release.py freed the input coordinates by zeroing their gate, and then measured the
Jacobian of (A, B) in them as identically zero.  That is correct and it is the point:
the coordinates that enter the addition are x12186, x16742, x14853, x24908, which
only EQUAL the released constants because the advice DAG puts them there.  Moving
x6418 does nothing until the DAG is re-solved -- `x6418 -> x1308 -> x14853`.

So the map to model is the composite

    F(released coordinates)  =  (A, B)   after re-running the advice sweep

which is still polynomial (§122), so fpoly interpolates it exactly.  Scan one
released coordinate for the roots of A, then the other for the roots of B, and
combine.  Everything measured, nothing linearised.

Usage: closer.py [gate] [state.json] [K]
"""
import os, sys, time
HERE = '/home/user/integer_solver/solve_lab/s10'
LAB = '/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
import fpoly as UP
import suppfree
from intad import jacZ
P = ad.P
gate = sys.argv[1] if len(sys.argv) > 1 else '24601'
src = sys.argv[2] if len(sys.argv) > 2 else 'PIN_39013.json'
K = int(sys.argv[3]) if len(sys.argv) > 3 else 8
GATES = {'2081': [2081], '24601': [24601], 'both': [2081, 24601]}[gate]
FREEDBY = {2081: [(6418, 'x2'), (12553, 'y2')], 24601: [(22152, 'x1'), (33462, 'y1')]}
KN = [t for g in GATES for t, _ in FREEDBY[g]]

v0 = L.load(src if os.path.isabs(src) else os.path.join(HERE, src))
ad.fwd(v0, rounds=6)
for g in GATES:
    v0[g] = 0
ad.fwd(v0, rounds=6)
print('%s, gate(s) %s released; knobs %s' % (src, GATES, KN), flush=True)

# the advice DAG, hard-coded from advgraph's output (topological order)
CONST = {6418: 20302955751113177691132960011219991444785130617995423281601414462835238472546,
         12553: 4531249068709477613185164105669741036354237152756954144434674493737552368539,
         22152: 82007976112976807461901870199198737303514020147647909878034348606308756230357,
         33462: 37841415183514949237467304684128824427406379377151921996714091976892367869714}
TWO = [(8778, 16144), (14623, 27522), (14853, 1308), (16742, 19083),
       (22649, 29524), (24548, 25442), (31339, 6858)]
ORDER_T = [6418, 33462, 8778, 12553, 24548, 14623, 14853, 16742, 22152, 22649, 31339]


def advsolve(v, held=()):
    v = list(v)
    for t, y in TWO:
        pass
    for _ in range(3):
        for t in ORDER_T:
            if t in held:
                continue
            tgt = None
            for a, y in TWO:
                if a == t:
                    tgt = v[y] % P
            if tgt is None:
                continue
            v[t] = (v[t] // P) * P + tgt
            ad.fwd(v, rounds=6)
    return v


def F(vals):
    v = list(v0)
    for u, x in vals.items():
        v[u] = x
    ad.fwd(v, rounds=6)
    v = advsolve(v, held=set(vals))
    return v, (v[35389] % P, v[6671] % P)


v, ab = F({})
print('base after advice re-solve: A = %d\n                            B = %d' % ab,
      flush=True)
print('score %d' % (L.NEQ - len(L.failing_eqs(L.all_atom_values(v)))), flush=True)

t0 = time.time()
for u in KN:
    xs = list(range(K + 3))
    ys = []
    for d in xs:
        _, ab2 = F({u: v0[u] + d})
        ys.append(ab2)
    fa = UP.interp(xs[:K + 1], [y[0] for y in ys])
    fb = UP.interp(xs[:K + 1], [y[1] for y in ys])
    oka = all(sum(c * pow(x, e, P) for e, c in enumerate(fa)) % P == ys[x][0]
              for x in xs[K + 1:])
    okb = all(sum(c * pow(x, e, P) for e, c in enumerate(fb)) % P == ys[x][1]
              for x in xs[K + 1:])
    print('\nknob x%d: deg(A) = %d exact=%s ; deg(B) = %d exact=%s  (%.0fs)'
          % (u, len(fa) - 1, oka, len(fb) - 1, okb, time.time() - t0), flush=True)
    for nm, f, ok in (('A', fa, oka), ('B', fb, okb)):
        if not ok or len(f) <= 1:
            continue
        rs = UP.roots(f)
        print('   %s has %d root(s) in F_p along this knob' % (nm, len(rs)), flush=True)
        for r in rs[:4]:
            w, ab3 = F({u: v0[u] + r})
            aw = L.all_atom_values(w)
            s = L.NEQ - len(L.failing_eqs(aw))
            nz = [a for a in range(L.NA) if a not in L.atom_out and aw[a]]
            print('      root -> A=%s B=%s  score %d  checks %d'
                  % (str(ab3[0])[:14], str(ab3[1])[:14], s, len(nz)), flush=True)
            T.save(w, os.path.join(HERE, 'CL_%d_x%d_%s.json' % (s, u, nm)))
    gcdf = UP.pgcd(fa, fb)
    if gcdf and len(gcdf) > 1:
        cr = UP.roots(gcdf)
        print('   COMMON roots of A and B along x%d: %d' % (u, len(cr)), flush=True)
        for r in cr[:4]:
            w, ab3 = F({u: v0[u] + r})
            aw = L.all_atom_values(w)
            s = L.NEQ - len(L.failing_eqs(aw))
            print('      *** common root -> score %d' % s, flush=True)
            T.save(w, os.path.join(HERE, 'CLC_%d.json' % s))
print('\ndone (%.0fs)' % (time.time() - t0))
