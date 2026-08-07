#!/usr/bin/env python3
"""K32: the degeneracy negative, redone properly for EVERY stage.

WHAT WAS WRONG BEFORE.  My interior-stage argument said "a stage with n < 256 leaves has
|x-y| < 2^n <= 2^255 < N, so x = y".  That is false: a stage's exponent set is an arbitrary
subset of {0..255}, not an initial segment, so a stage owning exponent 255 has x up to 2^255
and beyond regardless of how few leaves it has.  Only the ROOT case was argued correctly.
Agent P's challenge is what surfaced this.

THE MODULUS QUESTION, ANSWERED.  The modulus governing equality of the two coordinate pairs
is N (the order of the chain base, measured in k21 by composing the chain itself; N is prime
and the base is not the identity, so the order is exactly N).  N is NOT greater than 2^256 --
it is just below it.  The condition the carry walk actually needs is weaker:

    every stage has  |x - y| <= sum_{e in J1 u J2} 2^e  <  2^256  <  2N,

so the only multiples of N in range are  k = 0, +1, -1.   2N - 2^256 > 0 is checked below.
k = 0 needs x = y on disjoint bit supports, i.e. x = y = 0, excluded since both halves must be
live.  So k = +-1 exhausts the cases and a two-direction carry walk IS complete.

REDUCTION THAT MAKES THIS CHEAP AND EXACT.  |x-y| = N > 2^255 requires one side's exponent set
J to satisfy sum_{e in J} 2^e >= N > 2^255.  But sum_{e <= 254} 2^e = 2^255 - 1 < N.  So J must
contain exponent 255.  Only stages whose subtree contains that one leaf can possibly be
degenerate -- the ancestors of it, a chain of at most ~8 stages.  Those are enumerated here by
set inclusion, no parent pointers needed.

The per-bit rule now has THREE cases, because positions outside J1 u J2 are available to
neither side and can therefore FAIL mid-walk, not only at the final carry."""
import sys, os, json, collections
K = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, K)
F = '/home/user/integer_solver/solve_lab/agentF_work'
sys.path.insert(0, F)
import mux as MUX
import re

N = int(json.load(open(K + '/order.json'))['N'][0])
assert 2 * N > 2 ** 256 > N, 'the k=+-1-only argument needs 2N > 2**256'
print('N > 2^255 :', N > 2 ** 255, '   2N > 2^256 :', 2 * N > 2 ** 256, '  (so k in {0,+1,-1} only)')
print('sum_{e<=254} 2^e = 2^255-1 < N :', 2 ** 255 - 1 < N, ' => a degenerate side must own exponent 255')

ch = json.load(open(K + '/chain.json'))
sel2exp = {ch['sel'][str(i)]: ch['exp'][str(i)] for i in range(256)}
exp2sel = {e: s for s, e in sel2exp.items()}

# ---- exponent support of every decodable wire (same descent as K20) -------------------
import fold as FD
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

TOP = exp2sel[255]
sets = set(s for s in memo.values() if TOP in s and len(s) > 1)
sets.add(frozenset(exp2sel[e] for e in range(256)))          # the root, known exactly
chain = sorted(sets, key=len)
print('\nnested supports containing exponent 255 (ancestor chain), sizes:', [len(s) for s in chain])

MISSING = set(exp2sel.values()) - set().union(*chain) if chain else set()
print('leaves never assigned by the descent:', sorted(sel2exp[s] for s in MISSING))


def walk(m, plus, minus, allowed):
    """x (bits in `plus`) - y (bits in `minus`) == m, positions outside `allowed` forced 0 in
    BOTH.  Deterministic; returns None the moment a bit has no legal choice."""
    x = y = c = 0
    for i in range(256):
        mi = (m >> i) & 1
        if i in plus:
            s = mi + c; x |= (s & 1) << i; c = s >> 1
        elif i in minus:
            s = mi + c
            if s == 0: yi, c = 0, 0
            elif s == 1: yi, c = 1, 1
            else: yi, c = 0, 1
            y |= yi << i
        else:                       # neither side may use this bit
            s = mi + c
            if s == 1: return None  # x_i must be 0 but the sum is odd -> dead
            c = s >> 1
    return None if c else (x, y)


print('\n%-6s %-8s %-8s  x-y=N   y-x=N' % ('level', '|J1|', '|J2|'))
anyok = []
for k in range(len(chain)):
    S1 = chain[k]
    S2 = chain[k + 1] if k + 1 < len(chain) else None
    if S2 is None: break
    J1 = set(sel2exp[s] for s in S1)
    J2 = set(sel2exp[s] for s in (S2 - S1))
    variants = [(J1, J2)]
    for s in (MISSING & S2):        # conservative: try the unassigned leaf on either side
        e = sel2exp[s]
        variants += [(J1 | {e}, J2), (J1, J2 | {e})]
    res = []
    for j1, j2 in variants:
        al = j1 | j2
        res.append(walk(N, j1, j2, al) is not None or walk(N, j2, j1, al) is not None)
    ok = any(res)
    print('%-6d %-8d %-8d  %s' % (k, len(J1), len(J2), 'DEGENERACY POSSIBLE' if ok else 'no'))
    if ok: anyok.append((k, len(J1), len(J2)))

print('\nstages admitting a degeneracy:', anyok or 'NONE')
print('stages that cannot own exponent 255 are excluded by the size bound above, exactly.')
