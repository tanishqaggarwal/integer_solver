#!/usr/bin/env python3
"""K36: tighten the leaf-support recovery, then redo the disjoint-pair degeneracy test.

WHY.  K33's tree-free test found pairs admitting x-y=N and so FAILED to close the negative.
But K33's supports are over-approximations: for a gated mux term (xA*xB) my descent did not
know which operand was the selector and which was the value, so it unioned BOTH cones.  A
quadrant gate's selector is a liveness bit of the SIBLING subtree, so sibling leaves leak in.
That is why supports of size 252/253 appeared at all -- they cannot be children of a root that
splits 178/78.

FIX.  Identify liveness/boolean variables by fixpoint (free booleans, plus anything defined as
a product, sum or complement of them), and in a gated term take only the operand that is NOT
one of those.  Sanity gate: the root must come back 178/78, or the recovery is still wrong and
the test below is not to be believed."""
import sys, os, json, re, collections
K = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, K)
F = '/home/user/integer_solver/solve_lab/agentF_work'
sys.path.insert(0, F)
import mux as MUX
import fold as FD
from parse import node_str
from circ2 import vars_of

N = int(json.load(open(K + '/order.json'))['N'][0])
ch = json.load(open(K + '/chain.json'))
sel2exp = {ch['sel'][str(i)]: ch['exp'][str(i)] for i in range(256)}
exp2sel = {e: s for s, e in sel2exp.items()}
vc = json.load(open(K + '/varclass2.json'))

E = MUX.E
defnode = {}
for a in E.order:
    c = E.cls[a]; defnode[c[1]] = c[2]

# ---- LIVE = boolean / liveness variables, by fixpoint --------------------------------
LIVE = set(vc['bools'])
changed = True
while changed:
    changed = False
    for w, n in defnode.items():
        if w in LIVE: continue
        vs = vars_of(n) - {w}
        if vs and vs <= LIVE:
            s = node_str(n)
            if re.fullmatch(r'\(?x\d+\)?', s) or '*' in s or s.startswith('(1-') or '+' in s:
                LIVE.add(w); changed = True
print('liveness/boolean variables (fixpoint):', len(LIVE))

D = FD.points()
leafsel = {}
for l in D['leaves']:
    leafsel[l['wx']] = l['sel']; leafsel[l['wy']] = l['sel']
gatedpat = re.compile(r'^\(x(\d+)\*x(\d+)\)$')
memo = {}
AMBIG = [0]


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
                if not m: continue
                a, b = int(m.group(1)), int(m.group(2))
                va = [u for u in (a, b) if u not in LIVE]
                if len(va) != 1: AMBIG[0] += 1; va = [a, b]
                for u in va: out |= support(u, depth + 1)
            elif kind == 'free':
                out |= support(t, depth + 1)
    memo[w] = frozenset(out)
    return memo[w]


cands = set()
for a in E.res:
    for u in re.findall(r'x(\d+)', a): cands.add(int(u))
for w in sorted(cands):
    try: support(w)
    except Exception: pass
print('gated terms where the selector could not be told from the value:', AMBIG[0])

rootA = support(12186) | support(16742)
rootB = support(14853) | support(24908)
print('SANITY GATE  root A support = %d   root B support = %d   disjoint=%s   union=%d'
      % (len(rootA), len(rootB), not (rootA & rootB), len(rootA | rootB)))
OK = (len(rootA), len(rootB)) in ((178, 78), (177, 78), (178, 77)) and not (rootA & rootB)
print('recovery credible:', OK)

sets = set(s for s in memo.values() if s)
extra = set()
for a in sets:
    for b in sets:
        if a < b: extra.add(b - a)
sets |= {s for s in extra if s}
sets = sorted(sets, key=len)
expof = {s: frozenset(sel2exp[u] for u in s) for s in sets}
maskval = {s: sum(1 << e for e in expof[s]) for s in sets}
print('distinct supports:', len(sets))
plusside = [s for s in sets if maskval[s] >= N]
print('supports able to serve as the "+" side (maskval >= N):', len(plusside))


def walk(m, plus, minus):
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
        else:
            s = mi + c
            if s == 1: return None
            c = s >> 1
    return None if c else (x, y)


hits = []
tested = 0
for a in plusside:
    ja = expof[a]
    for b in sets:
        if a is b or (a & b): continue
        tested += 1
        if walk(N, ja, expof[b]) is not None:
            hits.append((len(ja), len(expof[b])))
print('disjoint pairs tested:', tested)
print('pairs admitting x - y = N :', len(hits), sorted(set(hits))[:20])
print('\nCONCLUSION:', 'no disjoint pair of recovered supports admits a degeneracy'
      if not hits else 'NOT CLOSED -- some pair admits one; must check if it is a real stage')
