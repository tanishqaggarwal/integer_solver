"""Corrected stage B: rank placements by INTEGER REACHABILITY of the row targets.

Differences from agentH_work/stageB.py:
  1. solvability of a row subset is tested COMPLETELY (zsolve, HNF over Z) instead of by
     the pivot-only particular solution, which produced false negatives;
  2. the maximum integrally-zeroable row count is found by downward-closed DFS instead of
     enumerating subsets of size <= min(n, rank, 8), so it is exact and not rank-capped;
  3. state construction is incremental, so a placement costs ~0.1 s rather than ~1 s.

Score model (identical to H's, so numbers are comparable):
    failing = |R| - optimum + outside ,  score = 39033 - failing
"""
import ev, json, time, sys, os
import frameB as FB
from frameB import Frame, State
from collections import defaultdict
from fractions import Fraction
import ast, re
import zsolve

VAR_RE = re.compile(r'x_(\d+)')
HERE = os.path.dirname(os.path.abspath(__file__))

atom_eqs = defaultdict(set)
for i, (m, sq, tl) in enumerate(ev.eq_terms):
    for c, a in tl:
        atom_eqs[a].add(i)

POOL = json.load(open(os.path.join(HERE, 'pool.json')))
WIT = [642, 28730, 29854, 31864]
ORDER = [v for v in ev.F['order'] if v in set(POOL)]
DEFRHS = {}
for v in POOL:
    a = ev.F['definer'][v]
    t = ast.parse(ev.atom_src[a], mode='eval').body
    DEFRHS[v] = compile(VAR_RE.sub(r'v[\1]', ast.unparse(t.right)), '<d>', 'eval')

_t0 = time.time()
fr = Frame(POOL)
FREE = set(fr.free)
FR0 = fr.free
W = json.load(open(os.path.join(HERE, '..', 'best', 'new_instance_partial_39026.json')))
wv = [0] * 38748
for k, val in W.items():
    wv[int(k[2:]) if k.startswith('x_') else int(k)] = int(val)
BASEFV = {u: wv[u] for u in fr.free if wv[u] != 0}


_BASE = [None]


def make(D, passes=1):
    """State with pool set D detached (holding witness values) and every other pool var
    re-attached to its gate value.  ORDER is topological, so one pass suffices; `passes`
    is kept for cross-checking against H's 3-pass version."""
    if _BASE[0] is None:
        _BASE[0] = State(fr, BASEFV)
    st = _BASE[0].clone()
    Ds = set(D)
    for _ in range(passes):
        for v in ORDER:
            if v in Ds:
                continue
            st.set_free({v: eval(DEFRHS[v], st.ns)})
    return st


def _bits(x):
    o = []
    while x:
        q = x & -x
        o.append(q.bit_length() - 1)
        x ^= q
    return o


def inner(st, e):
    m, sq, tl = ev.eq_terms[e]
    t = 0
    for c, a in tl:
        x = st.av.get(a)
        if x:
            t += c * x
    return t


def rank_q(rows, n):
    rr = [[Fraction(x) for x in r] for r in rows]
    piv = 0
    for c in range(n):
        k = None
        for i in range(piv, len(rr)):
            if rr[i][c] != 0:
                k = i
                break
        if k is None:
            continue
        rr[piv], rr[k] = rr[k], rr[piv]
        pv = rr[piv][c]
        for i in range(len(rr)):
            if i != piv and rr[i][c] != 0:
                f = rr[i][c] / pv
                rr[i] = [rr[i][j] - f * rr[piv][j] for j in range(n)]
        piv += 1
        if piv == len(rr):
            break
    return piv


def build(st):
    """Region R, cancel-set S, zero-collateral knobs, target b, knob-response matrix M."""
    NZ = set(st.nz())
    if not NZ:
        return None
    R = set()
    for q in NZ:
        R |= atom_eqs[q]
    Rl = sorted(R)
    S = sorted(q for q in set().union(*[set(z for c, z in ev.eq_terms[e][2]) for e in Rl])
               if atom_eqs[q] <= R)
    Sset = set(S)
    cands = set()
    for q in S:
        if q in fr.csup:
            cands.update(FR0[bb] for bb in _bits(fr.csup[q]))
    knobs = []
    for Y in sorted(cands):
        if Y not in FREE:
            continue
        h = st.clone().set_free({Y: st.fv.get(Y, 0) + 1})
        d = {q: h.av[q] - st.av[q] for q in h.av if h.av[q] != st.av[q]}
        if d and all(q in Sset for q in d):
            knobs.append(Y)
    b = [inner(st, e) for e in Rl]
    cols = []
    for Y in knobs:
        h = st.clone().set_free({Y: st.fv.get(Y, 0) + 1})
        cols.append([inner(h, e) - b[i] for i, e in enumerate(Rl)])
    n = len(knobs)
    M = [[cols[j][i] for j in range(n)] for i in range(len(Rl))]
    outside = len([e for e in st.fails if e not in R])
    return dict(R=Rl, S=S, knobs=knobs, b=b, M=M, n=n, outside=outside)


def price(st, node_cap=200000, want_rank=True, want_lin=True):
    d = build(st)
    if d is None:
        return dict(R=0, S=0, deficit=0, knobs=0, rank=0, z0=0, opt=0, outside=0,
                    failing=len(st.fails), score=39033 - len(st.fails), exhaustive=True,
                    rows=[], lin=True)
    Rl, M, b, n = d['R'], d['M'], d['b'], d['n']
    z0 = sum(1 for x in b if x == 0)
    rk = rank_q(M, n) if n else 0
    lin = None
    if want_lin and n:
        g = st.clone()
        for Y in d['knobs']:
            g.set_free({Y: st.fv.get(Y, 0) + 1})
        lin = all(inner(g, e) == b[i] + sum(M[i][j] for j in range(n)) for i, e in enumerate(Rl))
    opt, rows, exh, nodes = zsolve.max_zero_rows(M, b, n, len(Rl), node_cap=node_cap)
    fail = len(Rl) - opt + d['outside']
    return dict(R=len(Rl), S=len(d['S']), deficit=len(Rl) - len(d['S']), knobs=n, rank=rk,
                z0=z0, opt=opt, outside=d['outside'], failing=fail, score=39033 - fail,
                exhaustive=exh, rows=rows, nodes=nodes, lin=lin, knoblist=d['knobs'])


def price_D(D, **kw):
    return price(make(list(D)), **kw)


if __name__ == '__main__':
    print('frame built %.1fs' % (time.time() - _t0), flush=True)
    t = time.time()
    r = price_D(WIT)
    print('CALIBRATION WITNESS %s' % WIT, flush=True)
    print('  |R|=%d |S|=%d deficit=%d knobs=%d rank=%d z0=%d OPT=%d outside=%d failing=%d score=%d '
          'exhaustive=%s nodes=%d lin=%s  (%.2fs)'
          % (r['R'], r['S'], r['deficit'], r['knobs'], r['rank'], r['z0'], r['opt'],
             r['outside'], r['failing'], r['score'], r['exhaustive'], r['nodes'], r['lin'],
             time.time() - t), flush=True)
    print('  reproduces 39,026:', r['score'] == 39026, flush=True)
    # cross-check against H's stageB on a couple of sets it reported
    for D in ([17499], [17499, 20492], [642, 28730, 31864]):
        t = time.time()
        r = price_D(D)
        print('D=%-24s |R|=%-3d knobs=%-2d rank=%-2d z0=%-2d OPT=%-2d outside=%d failing=%-2d '
              'score=%d exh=%s (%.2fs)'
              % (str(D), r['R'], r['knobs'], r['rank'], r['z0'], r['opt'], r['outside'],
                 r['failing'], r['score'], r['exhaustive'], time.time() - t), flush=True)
