#!/usr/bin/env python3
"""AUDIT T32b -- WHAT blocks the |S|=17 residue.

At |S|=17 the two-wire pass finds bivariate roots that clear the surviving condition by direct
recomputation (250 of 250 verified) yet NONE passes the global guard.  That is the shared-wire
simultaneity shape, in contrast with |S|=8 where the very first verified root cleared with zero
collateral.  This names the collateral: which atoms go nonzero when the clearing shift is applied,
and whether it is always the same ones."""
import os, sys, json, itertools, random, collections
T = '/home/user/integer_solver/solve_lab/agentT_work'
L = os.path.join(T, 'mirror', 'L')
sys.path.insert(0, T)
import t_close2w as C          # reuses the same engine, fit2, pairroots
E = C.E; SL = C.SL; SHIFT = C.SHIFT; p = C.p; NV = C.NV
relift = C.relift; vars_of = C.vars_of; atomvalvars = C.atomvalvars
influences = C.influences; nzcount = C.nzcount; factor = C.factor; crt_list = C.crt_list

vv = [0]*NV
for k, val in json.load(open(os.path.join(T, 'close_T17w.json'))).items():
    vv[int(k[2:])] = int(val)
relift(vv)
r = E.run(vv)
nz = [E.res[i] for i, x in enumerate(r) if x]
print('|S|=17 two-wire end state: %d nonzero atoms' % len(nz))
for a in nz:
    print('    ', a[:100])
BASE = nzcount(vv)
VV0 = vv[:]
ZERO0 = set(i for i, x in enumerate(E.run(vv)) if x == 0)
NZ0 = set(i for i, x in enumerate(E.run(vv)) if x)
print('baseline global nonzero: %d' % BASE)

TGT = ('x24468', 'x18956')
viol = [a for a in SL if r[E.residx[a]] != 0 and SL[a] and r[E.residx[a]] % abs(SL[a]) != 0
        and not any(t in a for t in TGT)]
print('non-target violated conditions: %d  %s' % (len(viol), [a[:60] for a in viol]))
A = viol[0]; IA = E.residx[A]; Cc = abs(SL[A])//p
FC = sorted(factor(Cc).items())
print('\ncondition %s  c=%d=%s' % (A[:70], Cc, '*'.join('%d^%d' % qe for qe in FC)))
WS = sorted(set(x for x in vars_of(E.atoms[A]) if x in SHIFT) |
            set(x for x in atomvalvars[A] if x in SHIFT))
WSl = [w for w in WS if influences(vv, A, w)]
print('candidate wires: %d admitted of %d  %s   -> %d pair(s)'
      % (len(WSl), len(WS), WSl, len(WSl)*(len(WSl)-1)//2))

rnd = random.Random(5)
for w, v in itertools.combinations(WSl, 2):
    cf = C.fit2(vv, IA, w, v)
    if cf is None:
        print('pair (x%d,x%d): not p-divisible on the grid' % (w, v)); continue
    ok = all(C.probe2(vv, IA, w, tw, v, tv)//p == C.peval2_exact(cf, tw, tv)
             for tw, tv in ((7, 11), (1234, 99), (55555, 4)))
    cross = any(cf[k][l] for k in range(1, len(cf)) for l in range(1, len(cf)))
    per = []
    for q, e in FC:
        rs, ex = C.pairroots(cf, q, e, 24)
        per.append((q**e, rs, ex))
    print('pair (x%d,x%d) fit-vs-recompute %s  cross=%s  roots/prime %s'
          % (w, v, 'OK' if ok else 'MISMATCH', cross, [len(rs) for _, rs, _ in per]))
    if any(not rs for _, rs, _ in per):
        print('   NO ROOT'); continue
    mods = [m for m, _, _ in per]; sets = [rs for _, rs, _ in per]
    pat = collections.Counter(); seen = set(); tried = cleared = 0
    best = None
    while tried < 400:
        pick = [rs[rnd.randrange(len(rs))] for rs in sets]
        tw = crt_list([(x, m) for (x, _), m in zip(pick, mods)])
        tv = crt_list([(y, m) for (_, y), m in zip(pick, mods)])
        if tw is None or tv is None or (tw, tv) in seen:
            tried += 1; continue
        seen.add((tw, tv)); tried += 1
        if (tw, tv) == (0, 0) or C.probe2(vv, IA, w, tw, v, tv) % (Cc*p) != 0:
            continue
        cleared += 1
        vv[w] += p*tw; vv[v] += p*tv
        n = nzcount(vv)
        rn = E.run(vv)
        new = tuple(sorted(i for i in ZERO0 if rn[i]))
        gone = tuple(sorted(i for i in NZ0 if not rn[i]))
        vv[:] = VV0
        pat[(new, gone, n)] += 1
        if best is None or n < best:
            best = n
    print('   tried %d, verified-clearing %d, best global %s (base %d)'
          % (tried, cleared, best, BASE))
    for (new, gone, n), cnt in pat.most_common(6):
        print('      x%-4d  global->%d   newly nonzero %d: %s' % (
            cnt, n, len(new), [E.res[i][:64] for i in new[:3]]))
        print('              cleared %d: %s' % (len(gone), [E.res[i][:64] for i in gone[:3]]))
