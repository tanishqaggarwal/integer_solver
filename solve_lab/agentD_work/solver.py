"""Agent D general repair engine.

Move generation for a nonzero check atom c:
  * probe each candidate free knob u for its exact slope ds = c(u+1) - c(u) and
    confirm affineness with a second probe;
  * single-knob exact solve  (ds | c0);
  * two-knob Bezout solve    (gcd(ds_i, ds_j) | c0)  -- this is what closes the
    p-quantised "advice" checks: knob = the advice variable (slope K), partner =
    its solo handle (slope +-p).
Every proposed move is applied and scored EXACTLY; nothing is trusted.
"""
import collections, random, time, sys, math, json, os
import dlib as L
import engine2 as E
import rad

P = L.P


def slopes(st, c, knobs, affine_check=True):
    """Return {u: slope} for knobs where atom c is affine in u."""
    c0 = st.av[c]
    out = {}
    for u in knobs:
        b = st.v[u]
        r = st.apply({u: b + 1})
        c1 = st.av[c]
        st.revert(r)
        if c1 == c0:
            continue
        if affine_check:
            r = st.apply({u: b + 2})
            c2 = st.av[c]
            st.revert(r)
            if c2 - c1 != c1 - c0:
                continue
        out[u] = c1 - c0
    return out


def gen_moves(st, c, knobs=None, pairs=True, maxpairs=4000):
    """Yield seed-dicts that (predicted) zero atom c."""
    c0 = st.av[c]
    if c0 == 0:
        return []
    if knobs is None:
        knobs = sorted(rad.free_knobs(c, st.v))
    sl = slopes(st, c, knobs)
    moves = []
    for u, s in sl.items():
        if c0 % s == 0:
            moves.append({u: st.v[u] - c0 // s})
    if pairs:
        items = sorted(sl.items(), key=lambda t: abs(t[1]))
        n = len(items)
        cnt = 0
        for i in range(n):
            ui, si = items[i]
            for j in range(i + 1, n):
                uj, sj = items[j]
                g = math.gcd(si, sj)
                if c0 % g:
                    continue
                # solve si*a + sj*b = -c0
                a0, b0 = _bez(si, sj)
                k = -c0 // g
                a, b = a0 * k, b0 * k
                # reduce a modulo sj/g to keep numbers small
                m = sj // g
                if m:
                    t = a // m
                    a -= t * m
                    b += t * (si // g)
                moves.append({ui: st.v[ui] + a, uj: st.v[uj] + b})
                cnt += 1
                if cnt >= maxpairs:
                    return moves
    return moves


def _bez(a, b):
    old_r, r = a, b
    old_s, s = 1, 0
    old_t, t = 0, 1
    while r:
        q = old_r // r
        old_r, r = r, old_r - q * r
        old_s, s = s, old_s - q * s
        old_t, t = t, old_t - q * t
    return old_s, old_t


def best_move(st, c, knobs=None, verbose=False, topn=1):
    res = []
    for mv in gen_moves(st, c, knobs):
        r = st.apply(mv)
        ok = st.av[c] == 0
        res.append((st.score, ok, mv))
        st.revert(r)
    res.sort(key=lambda t: (-t[0], not t[1]))
    if verbose:
        for s, ok, mv in res[:5]:
            print(f'    score {s} closed={ok} {list(mv)[:3]}')
    return res[:topn]


def repair_loop(st, maxiter=60, verbose=True, rng=None, noise=0.0):
    rng = rng or random.Random(0)
    for it in range(maxiter):
        nz = st.nz()
        if not nz:
            return True
        best = (st.score, None)
        order = sorted(nz, key=lambda a: -len(L.atom2eq.get(a, {})))
        cands = []
        for c in order:
            for s, ok, mv in best_move(st, c, topn=3):
                cands.append((s, ok, c, mv))
        if not cands:
            if verbose:
                print(f'  it{it}: no moves, score={st.score} nz={len(nz)}')
            return False
        cands.sort(key=lambda t: (-t[0], not t[1]))
        if noise > 0 and rng.random() < noise:
            pick = rng.choice(cands)
        else:
            pick = cands[0]
        if pick[0] <= st.score and noise == 0:
            if verbose:
                print(f'  it{it}: local optimum score={st.score} nz={len(nz)} best_cand={pick[0]}')
            return False
        st.apply(pick[3])
        if verbose:
            print(f'  it{it}: score={st.score} nz={len(st.nz())} (fixed a{pick[2]})')
    return False


if __name__ == '__main__':
    path = sys.argv[1] if len(sys.argv) > 1 else 'D_state1.json'
    st = E.St(L.load(path))
    print('start', st.score, st.nz())
    t0 = time.time()
    ok = repair_loop(st, maxiter=40)
    print('done', st.score, st.nz(), f'{time.time()-t0:.1f}s')
    L.save(st.v, path.replace('.json', '_rep.json'))
