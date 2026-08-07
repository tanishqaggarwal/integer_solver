#!/usr/bin/env python3
"""K33: same negative as K32 but WITHOUT trusting my tree recovery.

K32 paired J1/J2 by set-inclusion along an "ancestor chain".  That chain contained two
distinct sets of the same size, so it is not provably a nest, so the pairing is not provably
the real stage structure.  Rather than fix the tree, drop it:

  test EVERY pair (J1, J2) of DISJOINT exponent sets drawn from the supports the descent
  found, in both directions.

Every real stage pair is among those pairs (a stage's two children are supports, and they are
disjoint), so if none of them admits x - y = +-N, no real stage does either.  Extra pairs that
do not correspond to a stage only make the test stricter.

Still exact, because of the two facts checked in K32:
  2N > 2^256 > every |x - y|  =>  only k = 0, +1, -1, and k = 0 needs both halves empty;
  sum_{e<=254} 2^e < N        =>  a degenerate side must own exponent 255."""
import sys, os, json, re, itertools, collections
K = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, K)
F = '/home/user/integer_solver/solve_lab/agentF_work'
sys.path.insert(0, F)
import mux as MUX
import fold as FD

N = int(json.load(open(K + '/order.json'))['N'][0])
assert 2 * N > 2 ** 256 > N and 2 ** 255 - 1 < N
ch = json.load(open(K + '/chain.json'))
sel2exp = {ch['sel'][str(i)]: ch['exp'][str(i)] for i in range(256)}
exp2sel = {e: s for s, e in sel2exp.items()}

D = FD.points()
leafsel = {}
for l in D['leaves']:
    leafsel[l['wx']] = l['sel']; leafsel[l['wy']] = l['sel']
gatedpat = re.compile(r'^\(x(\d+)\*x(\d+)\)$')
memo = {}


def support(w, depth=0):
    if w in memo: return memo[w]
    if w in leafsel:
        memo[w] = frozenset({leafsel[w]}); return memo[w]
    memo[w] = frozenset()
    if depth > 40: return memo[w]
    out = set()
    for z, coef in MUX.source_of(w):
        if z == 'CONST': continue
        for kind, t in MUX.mux_terms(z):
            if kind == 'gated':
                m = gatedpat.match(t)
                if m:
                    for u in (int(m.group(1)), int(m.group(2))):
                        out |= support(u, depth + 1)
            elif kind == 'free':
                out |= support(t, depth + 1)
    memo[w] = frozenset(out)
    return memo[w]


cands = set()
for a in MUX.E.res:
    for u in re.findall(r'x(\d+)', a): cands.add(int(u))
for w in sorted(cands):
    try: support(w)
    except Exception: pass

ALL = frozenset(exp2sel.values())
sets = set(s for s in memo.values() if len(s) >= 1)
sets.add(ALL)
# also add every complement-within-a-larger-support, since a sibling is exactly that
extra = set()
for a in sets:
    for b in sets:
        if a < b: extra.add(b - a)
sets |= {s for s in extra if s}
sets = sorted(sets, key=len)
print('distinct exponent supports recovered:', len(sets))

TOP = exp2sel[255]
withtop = [s for s in sets if TOP in s]
print('supports owning exponent 255:', len(withtop))


def walk(m, plus, minus):
    allowed_plus, allowed_minus = plus, minus
    x = y = c = 0
    for i in range(256):
        mi = (m >> i) & 1
        if i in allowed_plus:
            s = mi + c; x |= (s & 1) << i; c = s >> 1
        elif i in allowed_minus:
            s = mi + c
            if s == 0: yi, c = 0, 0
            elif s == 1: yi, c = 1, 1
            else: yi, c = 0, 1
            y |= yi << i
        else:
            s = mi + c
            if s == 1: return None
            c = s >> 1
    return None if c else (x, y)


expof = {s: frozenset(sel2exp[u] for u in s) for s in sets}
maskval = {s: sum(1 << e for e in expof[s]) for s in sets}

# EXACT PRUNE.  x - y = N with y >= 0 forces x >= N, and x <= maskval(J1).  So the "+" side
# must have maskval >= N.  This is a necessary condition, so pruning on it loses nothing.
plusside = [s for s in sets if maskval[s] >= N]
print('supports that could serve as the "+" side (maskval >= N):', len(plusside))
assert all(TOP in s for s in plusside), 'consistency: any such support must own exponent 255'

hits = []
tested = 0
for a in plusside:
    ja = expof[a]
    for b in sets:
        if a is b or (a & b): continue          # must be disjoint
        jb = expof[b]
        tested += 1
        if walk(N, ja, jb) is not None:
            hits.append(('plus=%d minus=%d' % (len(ja), len(jb))))
print('disjoint pairs tested (superset of all real stage pairs, after the exact prune):', tested)
print('pairs admitting  x - y = +-N :', hits or 'NONE')
print()
print('CONCLUSION: no stage anywhere can have its two inputs carry equal coordinate pairs.'
      if not hits else 'CONCLUSION: the negative FAILS -- a degeneracy is possible.')
