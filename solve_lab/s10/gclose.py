"""S11 step 84: the general gadget-closing loop.

build7.py drove all seven residual atoms of the deliverable to EXACTLY ZERO by
construction -- the first time that has happened -- leaving four nonzero checks
instance-wide and no broken gate at all:

    a7930  = 9367949 *(x24548 - x25442) - x7927     x7927  = p*handle,  x24548 FREE
    a29539 = 12846437*(x14853 - x1308 ) - x29967    x29967 = p*x30163,  x14853 FREE
    a40826, a41512   bundle checks that CONTAIN those two residuals

So the whole instance now rests on two congruences on two free 296-bit advice
inputs:  x24548 ≡ x25442 (mod p)  and  x14853 ≡ x1308 (mod p).

This module closes such a gadget generically and iterates:

  * find, for a nonzero check c, a free input u whose mod-p Jacobian entry is
    nonzero -- that is an OPERAND, and jumping its residue kills c mod p;
  * find a free input h whose mod-p entry is zero but whose exact integer entry
    (intad.jacZ, the only thing in this lab that can see handles) is a multiple of
    p -- that is the HANDLE, and it absorbs the quotient over Z;
  * apply both, forward-evaluate, and record what the residual became.

Each round is exact.  Whether the loop converges, cycles, or wanders is the
measurement.

Usage: gclose.py [state.json] [ROUNDS]
"""
import os, sys, time, json
HERE = '/home/user/integer_solver/solve_lab/s10'
LAB = '/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
from fwdad import jac_column
from intad import jacZ
import suppfree
P = ad.P
src = sys.argv[1] if len(sys.argv) > 1 else 'B7_39004.json'
ROUNDS = int(sys.argv[2]) if len(sys.argv) > 2 else 12
v = L.load(src if os.path.isabs(src) else os.path.join(HERE, src))
ad.fwd(v, rounds=6)
_, freelist, SVS = suppfree.build(v, modp=None)


def supp(c, v):
    m = suppfree.atom_supp(c, v, SVS, modp=None)
    return [freelist[i] for i in range(len(freelist)) if (m >> i) & 1]


def state(v):
    av = L.all_atom_values(v)
    return (L.NEQ - len(L.failing_eqs(av)), av,
            [a for a in range(L.NA) if a not in L.atom_out and av[a]],
            [a for a in L.atom_out if av[a]])


s, av, NZ, GB = state(v)
print('%s: score %d; nonzero checks %s; broken gates %d'
      % (src, s, NZ, len(GB)), flush=True)
best, bestv = s, list(v)
seen = {}
for rnd in range(ROUNDS):
    if not NZ:
        print('*** NO NONZERO CHECK REMAINS'); break
    vm = [x % P for x in v]
    moved = False
    for c in NZ:
        sup = supp(c, v)
        ops = [(u, jac_column(u, v, vm, [c]).get(c, 0)) for u in sup]
        ops = [(u, g) for u, g in ops if g]
        handles = []
        for u in sup:
            if any(u == w for w, _ in ops):
                continue
            gz = jacZ(u, v, [c]).get(c, 0)
            if gz and gz % P == 0:
                handles.append((u, gz))
        if not ops:
            continue
        # jump the operand whose column touches the fewest other checks
        cand = sorted(ops, key=lambda t: len(jac_column(t[0], v, vm,
                      [a for a in range(L.NA) if a not in L.atom_out])))
        u, g = cand[0]
        d = (-av[c] % P) * pow(g, -1, P) % P
        w = list(v)
        w[u] = w[u] + d
        ad.fwd(w, rounds=6)
        aw = L.all_atom_values(w)
        if aw[c] % P:
            continue
        for h, gz in handles:
            if aw[c] % gz == 0:
                w[h] = w[h] - aw[c] // gz
                ad.fwd(w, rounds=6)
                aw = L.all_atom_values(w)
                break
        s2 = L.NEQ - len(L.failing_eqs(aw))
        nz2 = [a for a in range(L.NA) if a not in L.atom_out and aw[a]]
        print('  round %-2d close a%-6d via x%-6d (%d handles): score %d -> %d, '
              'checks %s' % (rnd, c, u, len(handles), s, s2, nz2), flush=True)
        v, av, NZ, s = w, aw, nz2, s2
        moved = True
        if s2 > best:
            best, bestv = s2, list(w)
            T.save(w, os.path.join(HERE, 'GC_%d.json' % s2))
            print('    *** saved GC_%d.json' % s2, flush=True)
        break
    if not moved:
        print('  round %d: no closable gadget' % rnd)
        break
    key = tuple(sorted(NZ))
    if key in seen:
        print('  CYCLE: residual set %s repeats (round %d, first seen %d)'
              % (list(key), rnd, seen[key]))
        break
    seen[key] = rnd
print('\nbest reached: %d' % best)
