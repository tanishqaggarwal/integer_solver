"""S11 step 112: is the cost of the two COARSE generators repairable?

The zero-cost lattice has nine generators and, measured exactly, no subset of more
than five of the twelve rows is integrally reachable inside it.  The obstruction is
the two coarse directions:

    a0 = a22229 moves only in multiples of 7376877   (x642, coupled to a6)
    a1 = a22230 moves only in multiples of p          (x9413)

and one of the twelve rows is `eq 29125 = {a22230: 1}` -- it demands a1 = 0 exactly,
which the lattice cannot reach because a1 is not ≡ 0 (mod p) at the witness.

The fine generators for those directions are x7068 (measured cost 13) and x28730
(measured cost 16), and with both the lattice is complete: alpha = 0 becomes
reachable, all twelve rows hold, and the score is 39,033 - 29 = 39,004 -- exactly
what build7 measured.  So the whole question is whether those 13 and 16 broken
equations can be REPAIRED by other free inputs.  If they can, the lattice is complete
at zero cost and every one of the twelve rows is reachable.

This module moves the coarse generator, lists precisely what breaks, and tries to
repair each broken equation with the machinery that exists: the integer lift for the
ones that are ≡ 0 (mod p), and a single-knob exact solve for the rest.

Usage: repaircost.py [7068|28730] [DELTA]
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
who = int(sys.argv[1]) if len(sys.argv) > 1 else 7068
DELTA = int(sys.argv[2]) if len(sys.argv) > 2 else 1
SEVEN = [22229, 22230, 35758, 35759, 35760, 35761, 35762]
base = L.load(os.path.join(LAB, 'best', 'new_instance_partial_39026.json'))
fwd(base)
E12 = set(e for a in SEVEN for e in L.atom2eq[a])
av0 = L.all_atom_values(base)
print('witness %d' % (L.NEQ - len(L.failing_eqs(av0))), flush=True)

v = list(base)
v[who] = v[who] + DELTA
fwd(v)
av = L.all_atom_values(v)
broke = sorted(set(L.failing_eqs(av)) - E12)
print('x%d += %d  ->  %d equations break outside the twelve: %s'
      % (who, DELTA, len(broke), broke), flush=True)


def comb(a_vals, e):
    s = 0
    for a, c in L.eq_atoms[e][2].items():
        if a_vals[a]:
            s += c * a_vals[a]
    return s


mp = [e for e in broke if comb(av, e) % P == 0]
print('   of those, combination ≡ 0 (mod p) -- liftable: %d %s'
      % (len(mp), mp), flush=True)
nzc = sorted({a for e in broke for a in L.eq_atoms[e][2] if av[a]})
print('   nonzero atoms in the broken equations: %s' % nzc[:20], flush=True)

_, freelist, SVS = suppfree.build(v, definer=definer, ORDER=ORDER, FREE=FREE,
                                  modp=None)
t0 = time.time()
cur = L.NEQ - len(L.failing_eqs(av))
print('\nscore after the move: %d; trying to repair each broken equation'
      % cur, flush=True)
for it in range(40):
    av = L.all_atom_values(v)
    broke = sorted(set(L.failing_eqs(av)) - E12)
    if not broke:
        print('   *** every collateral equation repaired', flush=True)
        break
    cur = L.NEQ - len(L.failing_eqs(av))
    best = (cur, None, None, None)
    for e in broke[:6]:
        S = comb(av, e)
        atoms = sorted(L.eq_atoms[e][2])
        m = 0
        for a in atoms:
            m |= suppfree.atom_supp(a, v, SVS, modp=None)
        cand = [freelist[i] for i in range(len(freelist)) if (m >> i) & 1]
        for u in cand[:120]:
            col = jacZ(u, v, atoms)
            g = 0
            for a, c in L.eq_atoms[e][2].items():
                if a in col:
                    g += c * col[a]
            if not g or S % g:
                continue
            w = list(v)
            w[u] = w[u] - S // g
            fwd(w)
            aw = L.all_atom_values(w)
            s2 = L.NEQ - len(L.failing_eqs(aw))
            if comb(aw, e) == 0 and s2 > best[0]:
                best = (s2, w, e, u)
    if best[1] is None:
        print('   no strictly-improving repair at score %d (%d still broken)'
              % (cur, len(broke)), flush=True)
        break
    cur, v, e, u = best
    print('   eq %d repaired by x%d -> score %d  (%.0fs)'
          % (e, u, cur, time.time() - t0), flush=True)
av = L.all_atom_values(v)
s = L.NEQ - len(L.failing_eqs(av))
broke = sorted(set(L.failing_eqs(av)) - E12)
print('\nFINAL score %d; %d collateral equations still broken: %s'
      % (s, len(broke), broke[:12]))
T.save(v, os.path.join(HERE, 'RC_%d_%d.json' % (who, s)))
print('saved RC_%d_%d.json' % (who, s))
