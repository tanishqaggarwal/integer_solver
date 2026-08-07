"""Corrected stage B over the OTHER two carrier classes: 20 cascade pins, 1,147 handles.

Agent H ran stage B only on top scorers under the (later withdrawn) rank>deficit criterion.
Here every carrier gets the exact integer-reachability number from zsolve.max_zero_rows.

usage: python3 carrierN.py pins
       python3 carrierN.py handles <shard> <nshard>
"""
import ev, fast, json, time, sys, os
from fast import St, csup
from chain import close_trace
from collections import defaultdict
from fractions import Fraction
import zsolve

HERE = os.path.dirname(os.path.abspath(__file__))
atom_eqs = defaultdict(set)
for i, (m, sq, tl) in enumerate(ev.eq_terms):
    for c, a in tl:
        atom_eqs[a].add(i)
BITS = json.load(open(os.path.join(HERE, 'bits.json')))
ALL = set(BITS['A'] + BITS['B'])
FREE = set(ev.F['free0'])
FR0 = ev.F['free0']


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


def price_state(out, node_cap=120000):
    NZ = set(out.nz())
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
        if q in csup:
            cands.update(FR0[bb] for bb in _bits(csup[q]))
    knobs = []
    for Y in sorted(cands):
        if Y not in FREE:
            continue
        h = out.clone().set_free({Y: out.fv.get(Y, 0) + 1})
        d = {q: h.av[q] - out.av[q] for q in h.av if h.av[q] != out.av[q]}
        if d and all(q in Sset for q in d):
            knobs.append(Y)
    b = [inner(out, e) for e in Rl]
    cols = []
    for Y in knobs:
        h = out.clone().set_free({Y: out.fv.get(Y, 0) + 1})
        cols.append([inner(h, e) - b[i] for i, e in enumerate(Rl)])
    n = len(knobs)
    M = [[cols[j][i] for j in range(n)] for i in range(len(Rl))]
    rk = rank_q(M, n) if n else 0
    z0 = sum(1 for x in b if x == 0)
    outside = len([e for e in out.fails if e not in R])
    opt, rows, exh, tests = zsolve.max_zero_rows(M, b, n, len(Rl), node_cap=node_cap)
    fail = len(Rl) - opt + outside
    return dict(R=len(Rl), S=len(S), deficit=len(Rl) - len(S), knobs=n, rank=rk, z0=z0,
                opt=opt, outside=outside, failing=fail, score=39033 - fail, exh=exh,
                tests=tests)


if __name__ == '__main__':
    what = sys.argv[1]
    st0 = St({})
    bsel = BITS['A'][0]
    if what == 'pins':
        CHAIN = [a for a, X in json.load(open(os.path.join(HERE, 'chain.json')))]
        out = open(os.path.join(HERE, 'runs', 'pins.jsonl'), 'w')
        t0 = time.time()
        for P in CHAIN:
            st = st0.clone().set_free({bsel: 1})
            o, ok, tr, fr = close_trace(st, frozen=set(ALL), skip=frozenset([P]))
            r = price_state(o)
            if r is None:
                continue
            r['pin'] = P
            out.write(json.dumps(r) + '\n')
            out.flush()
            print('a%-7d |R|=%-3d knobs=%-2d rank=%-2d z0=%-2d OPT=%-3d outside=%d failing=%-3d '
                  'score=%d exh=%s tests=%d  %.0fs'
                  % (P, r['R'], r['knobs'], r['rank'], r['z0'], r['opt'], r['outside'],
                     r['failing'], r['score'], r['exh'], r['tests'], time.time() - t0), flush=True)
        out.close()
    else:
        shard = int(sys.argv[2]); nsh = int(sys.argv[3])
        H = [(int(X), int(a)) for X, a, m in json.load(open(os.path.join(HERE, 'handles.json')))]
        NAMED = [(7497, 688), (11436, 1618), (22820, 30982), (14393, 30980)]
        TARGETS = NAMED + H
        TARGETS = [t for i, t in enumerate(TARGETS) if i % nsh == shard]
        base, ok, tr, frz = close_trace(st0.clone().set_free({bsel: 1}), frozen=set(ALL))
        print('base closed state score', base.score(), 'targets', len(TARGETS), flush=True)
        path = os.path.join(HERE, 'runs', 'handles_%d.jsonl' % shard)
        done = set()
        if os.path.exists(path):
            for ln in open(path):
                try:
                    done.add(json.loads(ln)['handle'])
                except Exception:
                    pass
        out = open(path, 'a')
        t0 = time.time()
        bestsc = 0
        for i, (X, a) in enumerate(TARGETS):
            if X in done:
                continue
            try:
                g = base.clone().set_free({X: base.fv.get(X, 0) + 1})
                o, ok2, tr2, _ = close_trace(g, frozen=set(ALL) | {X})
                r = price_state(o)
            except Exception as e:
                out.write(json.dumps({'handle': X, 'err': str(e)[:200]}) + '\n')
                continue
            if r is None:
                continue
            r['handle'] = X
            r['carrier'] = a
            out.write(json.dumps(r) + '\n')
            if r['score'] > bestsc:
                bestsc = r['score']
            if r['score'] > 39026:
                print('*** BEATS 39026: handle x_%d score=%d ***' % (X, r['score']), flush=True)
            if (i + 1) % 25 == 0:
                out.flush()
                el = time.time() - t0
                print('%d/%d  %.2fs/carrier  eta %.1f min  best=%d'
                      % (i + 1, len(TARGETS), el / (i + 1), (len(TARGETS) - i - 1) * el / (i + 1) / 60,
                         bestsc), flush=True)
        out.flush()
        out.close()
        print('DONE shard %d  best=%d  %.1f min' % (shard, bestsc, (time.time() - t0) / 60), flush=True)
