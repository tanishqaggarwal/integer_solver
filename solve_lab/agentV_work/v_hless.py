#!/usr/bin/env python3
"""agent V -- V4.  STRUCTURAL CENSUS OF THE HANDLE-LESS ATOM POPULATION.

An atom's "handle" is a free variable in its cone that moves it by an exact multiple of p, so the
lift can absorb any residual into it.  `relift` is indexed by handle; `SL` (the slope table) only
has rows for atoms that HAVE one.  An atom with NO handle has nothing to absorb into: it must be
EXACTLY ZERO over Z, not merely divisible by something.  closeS4's machinery never sees these.

This file establishes WHAT THE POPULATION IS -- not whether a particular solver can reach it.

Reported:
  1. the exact split, and whether "handle-less" is a shape property, a position property, or a
     property of the guard that carries the atom;
  2. what free variables they actually depend on (this is the mechanism: `atomh` is computed over
     the atom's FREE-VARIABLE CLOSURE, so handle-less means no handle appears anywhere in the cone);
  3. how many of them are structurally CONSTANT (identically zero for every assignment) versus
     genuinely live;
  4. whether they carry SHIFT wires at all -- an atom no shift can move is not a solver gap.

Usage: python3 v_hless.py
"""
import os, sys, json, re, collections, random, time
import v_base as B

V = '/home/user/integer_solver/solve_lab/agentV_work'
E, SL, p, SHIFT, NV = B.E, B.SL, B.p, B.SHIFT, B.NV
g = B.GL
atomh = g['atomh']; handle = g['handle']; fa = g['fa']; vars_of = g['vars_of']
defrhs = g['defrhs']; M = g['M']

HL = [a for a in E.res if len(atomh[a]) == 0]
HD = [a for a in E.res if len(atomh[a]) == 1]


def shape(a):
    """coarse syntactic class of an atom, from its own text."""
    s = re.sub(r'x\d+', 'X', a)
    s = re.sub(r'\d{2,}', 'K', s)
    return s


def freeset(a):
    s = set()
    for u in vars_of(E.atoms[a]):
        s |= fa(u)
    return s


if __name__ == '__main__':
    t0 = time.time()
    print('ATOMS %d = %d handle-less + %d with exactly one handle (none has two)'
          % (len(E.res), len(HL), len(HD)), flush=True)

    # ---------------- 1. shape
    print('\n--- 1. SHAPE: is "no handle" visible in the atom text? ---', flush=True)
    shl = collections.Counter(shape(a) for a in HL)
    shd = collections.Counter(shape(a) for a in HD)
    allsh = set(shl) | set(shd)
    print('%-34s %8s %8s  %s' % ('shape', 'handle-', 'handle+', 'verdict'), flush=True)
    shared = 0
    for s in sorted(allsh, key=lambda s: -(shl[s]+shd[s]))[:16]:
        v = ('ONLY handle-less' if shd[s] == 0 else
             'ONLY handled' if shl[s] == 0 else 'BOTH')
        if shl[s] and shd[s]:
            shared += shl[s]
        print('%-34s %8d %8d  %s' % (s[:34], shl[s], shd[s], v), flush=True)
    both = sum(shl[s] for s in allsh if shl[s] and shd[s])
    print('handle-less atoms whose SHAPE also occurs among handled atoms: %d of %d (%.0f%%)'
          % (both, len(HL), 100.0*both/len(HL)), flush=True)

    # ---------------- 2. free-variable dependence -- the mechanism
    print('\n--- 2. FREE-VARIABLE CONE: the actual mechanism ---', flush=True)
    live = set(M['live']); dead = set(M['dead'])
    cats = collections.Counter()
    nfree = collections.Counter()
    for a in HL:
        fs = freeset(a)
        nfree[min(len(fs), 9)] += 1
        if not fs:
            cats['NO free variable at all (constant atom)'] += 1
        elif fs & handle:
            cats['(impossible) contains a handle'] += 1
        elif fs & live:
            cats['depends on leaf selectors'] += 1
        elif fs & SHIFT:
            cats['depends on a SHIFT wire but no leaf'] += 1
        else:
            cats['other free vars only'] += 1
    for k, v in cats.most_common():
        print('   %-46s %d' % (k, v), flush=True)
    print('   free-variable-count histogram (9+ pooled): %s' % dict(sorted(nfree.items())), flush=True)

    # ---------------- 3. structurally constant vs live
    print('\n--- 3. ARE THEY EVER NONZERO?  (random assignments, direct recomputation) ---',
          flush=True)
    rnd = random.Random(4242)
    everynz = collections.Counter()
    seen_nz = set()
    trials = []
    for k in range(6):
        n = [1, 2, 8, 17, 32, 64][k]
        S = rnd.sample(M['live'], n)
        vv = B.greedy_init(S)
        r = E.run(vv)
        nz = set(E.res[i] for i, x in enumerate(r) if x)
        hlnz = nz & set(HL)
        seen_nz |= hlnz
        trials.append((n, len(nz), len(hlnz)))
        print('   |S|=%-3d  global nonzero %-4d  of which HANDLE-LESS %d'
              % (n, len(nz), len(hlnz)), flush=True)
    print('   distinct handle-less atoms seen nonzero across those 6 states: %d of %d'
          % (len(seen_nz), len(HL)), flush=True)

    # ---------------- 4. reachability: do they carry a SHIFT wire at all?
    print('\n--- 4. CAN A SHIFT EVEN MOVE THEM? ---', flush=True)
    withshift = [a for a in HL if B.wires_of(a)]
    print('   handle-less atoms with at least one SHIFT wire in scope: %d of %d'
          % (len(withshift), len(HL)), flush=True)
    nwh = collections.Counter(min(len(B.wires_of(a)), 9) for a in HL)
    print('   SHIFT-wire-count histogram (9+ pooled): %s' % dict(sorted(nwh.items())), flush=True)

    json.dump({'n_atoms': len(E.res), 'n_handleless': len(HL), 'n_handled': len(HD),
               'shape_shared': both,
               'categories': dict(cats),
               'trials': trials,
               'ever_nonzero': sorted(seen_nz),
               'with_shift_wire': len(withshift)},
              open(os.path.join(V, 'v_hless.json'), 'w'), indent=1)
    print('\nwall %.1f s -> v_hless.json' % (time.time()-t0), flush=True)
