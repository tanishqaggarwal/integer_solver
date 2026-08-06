"""S11 step 106: the lift belongs at the EQUATION level, not the atom level.

ceil15 exposed a gap.  At PF_best_39015 the mod-p coset count says 16 equations fail,
the checker says 18: two equations have an atom COMBINATION that vanishes mod p but
not over Z.  Every lift in this lab has worked atom-by-atom -- find a check that is
≡ 0 (mod p) and absorb it with a handle -- and an equation whose combination vanishes
while its individual atoms do not is invisible to that.

The equation-level version is the right object, and it is the same arithmetic one
level up: for a failing equation e whose combination S_e is ≡ 0 (mod p),

    dS_e/du  =  sum_a c_a * (da/du)_Z          exact, via intad.jacZ
    u <- u - S_e / (dS_e/du)                   whenever that coefficient divides it

and S_e becomes exactly zero over Z without any atom of e having to be zero.  That is
precisely the mechanism the 39,026 deliverable exploits (§152) -- it just does it by
luck of search rather than by construction.

Usage: eqlift.py [state.json]
"""
import os, sys, time
HERE = '/home/user/integer_solver/solve_lab/s10'
LAB = '/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
from intad import jacZ, dpartZ
import suppfree
P = ad.P
src = sys.argv[1] if len(sys.argv) > 1 else 'PF_best_39015.json'
v = L.load(src if os.path.isabs(src) else os.path.join(HERE, src))
ad.fwd(v, rounds=6)


def comb(av, e):
    s = 0
    for a, c in L.eq_atoms[e][2].items():
        if av[a]:
            s += c * av[a]
    return s


def report(v, tag):
    av = L.all_atom_values(v)
    f = sorted(L.failing_eqs(av))
    mp = [e for e in f if comb(av, e) % P == 0]
    print('%-30s score %-6d failing %-3d liftable (≡0 mod p) %d %s'
          % (tag, L.NEQ - len(f), len(f), len(mp), mp[:8]), flush=True)
    return L.NEQ - len(f), av, mp


s0, av, LIFTABLE = report(v, 'start (%s)' % os.path.basename(src))
_, freelist, SVS = suppfree.build(v, modp=None)
t0 = time.time()
for rnd in range(30):
    av = L.all_atom_values(v)
    cur = L.NEQ - len(L.failing_eqs(av))
    todo = [e for e in sorted(L.failing_eqs(av)) if comb(av, e) % P == 0]
    if not todo:
        break
    # Scan EVERY (equation, knob) pair and keep the STRICTLY best move.  Accepting
    # anything that merely does not lose oscillates forever: fixing eq 7469 with the
    # handle x30317 breaks eq 7123 through the same handle, and back again.
    bestmove = (cur, None, None, None)
    for e in todo:
        S = comb(av, e)
        atoms = sorted(L.eq_atoms[e][2])
        m = 0
        for a in atoms:
            m |= suppfree.atom_supp(a, v, SVS, modp=None)
        cand = [freelist[i] for i in range(len(freelist)) if (m >> i) & 1]
        for u in cand:
            col = jacZ(u, v, atoms)
            g = 0
            for a, c in L.eq_atoms[e][2].items():
                if a in col:
                    g += c * col[a]
            if not g or S % g:
                continue
            w = list(v)
            w[u] = w[u] - S // g
            ad.fwd(w, rounds=6)
            aw = L.all_atom_values(w)
            s2 = L.NEQ - len(L.failing_eqs(aw))
            if comb(aw, e) == 0 and s2 > bestmove[0]:
                bestmove = (s2, w, e, u)
    moved = bestmove[1] is not None
    if moved:
        cur, v, e, u = bestmove
        print('   eq %d absorbed by x%d -> score %d' % (e, u, cur), flush=True)
    if not moved:
        print('   no equation-level lift available this round', flush=True)
        break
s, av, mp = report(v, 'after the equation lift')
if s > s0:
    T.save(v, os.path.join(HERE, 'EL_%d.json' % s))
    print('saved EL_%d.json  (%.0fs)' % (s, time.time() - t0))
