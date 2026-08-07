"""Agent B: scan every free input; measure (cost in equations, effect on the four
residual congruences).  Fast: forward-evaluate then incrementally rescore."""
import os, sys, json, pickle, collections, time
os.environ.setdefault('ORIENT', 'orient3.pkl')
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentB_work')
import beval as E

P = 115792089237316195423570985008687907853269984665640564039457584007908834671663
Q = 7376877 * P
NV = 38748
facs, atoms, eqs = E.facs, E.atoms, E.eqs

occ = collections.defaultdict(list)
for i, p in enumerate(facs):
    vs = set()
    for m in p: vs.update(m)
    for v in vs: occ[v].append(i)
fac_atoms = collections.defaultdict(list)
for i, a in enumerate(atoms):
    for f in set(a): fac_atoms[f].append(i)
atom_eqs = collections.defaultdict(list)
for e, (sc, L, k) in enumerate(eqs):
    for c, a in L: atom_eqs[a].append(e)

def full_state(val):
    fv = [E.fval(f, val) for f in range(len(facs))]
    av = []
    for a in atoms:
        t = 1
        for f in a:
            t *= fv[f]
            if t == 0: break
        av.append(t)
    eqv = [sum(c*av[a] for c, a in L) for sc, L, k in eqs]
    return fv, av, eqv

def incr(val, changed, fv, av, eqv):
    dirty_f = set()
    for v in changed: dirty_f.update(occ[v])
    fv = dict((f, fv[f]) for f in dirty_f) if False else fv
    newf = {}
    dirty_a = set()
    for f in dirty_f:
        newf[f] = E.fval(f, val); dirty_a.update(fac_atoms[f])
    def getf(f): return newf.get(f, fv[f])
    newa = {}; dirty_e = set()
    for a in dirty_a:
        t = 1
        for f in atoms[a]:
            t *= getf(f)
            if t == 0: break
        newa[a] = t; dirty_e.update(atom_eqs[a])
    def geta(a): return newa.get(a, av[a])
    bad = 0
    for e in dirty_e:
        s = sum(c*geta(a) for c, a in eqs[e][1])
        if s != 0 and eqv[e] == 0: bad += 1
        elif s == 0 and eqv[e] != 0: bad -= 1
    return bad, newf, newa, dirty_e

def main():
    base = '/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json'
    v0 = E.load(base)
    fv0, av0, eqv0 = full_state(v0)
    ok0 = sum(1 for x in eqv0 if x == 0)
    print("base %d" % ok0)
    r0 = (v0[28730] % P, v0[8731] % P, v0[9118] % P, (v0[7068]-v0[2099]) % Q)
    print("base residues nonzero:", [i for i, x in enumerate(r0) if x])
    freeval0 = {v: v0[v] for v in E.free}
    deltas = [1, P, 7376877, Q]
    res = []
    t0 = time.time()
    for n, v in enumerate(E.free):
        for d in deltas:
            fvl = dict(freeval0); fvl[v] = fvl[v] + d
            val, nd, nf2 = E.forward(fvl, default=v0)
            changed = [u for u in range(NV) if val[u] != v0[u]]
            if not changed: continue
            bad, _, _, _ = incr(val, changed, fv0, av0, eqv0)
            r = (val[28730] % P, val[8731] % P, val[9118] % P, (val[7068]-val[2099]) % Q)
            moved = tuple(int(r[i] != r0[i]) for i in range(4))
            if bad <= 0 or any(moved):
                res.append((v, d, bad, moved, len(changed), nd))
        if (n+1) % 400 == 0:
            print("  scanned %d/%d  %.1fs  hits=%d" % (n+1, len(E.free), time.time()-t0, len(res)), flush=True)
    print("total interesting knobs:", len(res))
    pickle.dump(res, open('/home/user/integer_solver/solve_lab/agentB_work/knobs.pkl','wb'), -1)
    zero = [r for r in res if r[2] <= 0]
    print("zero-or-negative-cost knobs:", len(zero))
    mv = [r for r in res if any(r[3])]
    print("knobs that move a residue:", len(mv))
    for r in sorted(mv, key=lambda t: t[2])[:40]:
        print("   x%-6d d=%s cost=%+d moves=%s nchanged=%d" % (r[0], 'p' if r[1] == P else ('Q' if r[1] == Q else r[1]), r[2], r[3], r[4]))

if __name__ == '__main__':
    main()
