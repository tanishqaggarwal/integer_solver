"""Analysis of how the sacrifice set S can be grown to unlock more confined atoms.

Key accounting (verified empirically):  failing = |S| - maxzero(S), and
maxzero(S) = D - C where D = rank of the achievable atom-value lattice and
C = number of independent mod-P congruences it carries.
So a growth step only pays if it unlocks MORE atoms than equations it costs.
"""
import pickle, collections, sys, itertools, time
import lib as L, model as MD

v0 = L.load(L.BEST24)
MD.BASEP = [MD.prim_val(a, v0) for a in range(L.NA)]
S13 = frozenset([2554, 6816, 8124, 8680, 9123, 9421, 12231, 12270, 12350, 14584, 18673, 22044, 29125])


def extra_map(S):
    """atom -> frozenset(equations outside S that mention it)"""
    out = {}
    for a, d in L.atom2eq.items():
        e = frozenset(d) - S
        out[a] = e
    return out


def report(S, label='', quiet=False):
    S = frozenset(S)
    A = MD.confined_atoms(S)
    mod = MD.build(S, v0, verbose=False)
    D = rank_of(mod)
    if not quiet:
        print(f'{label} |S|={len(S)} |A|={len(A)} knobs={len(mod["knobs"])} D={D}')
    return mod, A, D


def rank_of(mod):
    """Rational rank of the knob->atom map."""
    A = mod['A']
    knobs = mod['knobs']
    rows = [[mod['M'][a].get(x, 0) for x in knobs] for a in A]
    from fractions import Fraction
    m = [[Fraction(c) for c in r] for r in rows]
    r = 0
    nr, nc = len(m), (len(knobs) if knobs else 0)
    for c in range(nc):
        piv = None
        for i in range(r, nr):
            if m[i][c]:
                piv = i
                break
        if piv is None:
            continue
        m[r], m[piv] = m[piv], m[r]
        pv = m[r][c]
        for i in range(nr):
            if i != r and m[i][c]:
                f = m[i][c] / pv
                for j in range(c, nc):
                    m[i][j] -= f * m[r][j]
        r += 1
    return r


if __name__ == '__main__':
    ex = extra_map(S13)
    # which single equations unlock new atoms?
    unlock1 = collections.defaultdict(list)
    for a, e in ex.items():
        if len(e) == 1:
            unlock1[next(iter(e))].append(a)
    print('single equations that unlock >=1 new atom:', len(unlock1))
    for e, ats in sorted(unlock1.items(), key=lambda kv: -len(kv[1]))[:25]:
        print(f'   +eq {e}: unlocks atoms {ats}')
    cnt = collections.Counter(len(e) for e in ex.values())
    print('atom |extra| histogram (<=6):', {k: cnt[k] for k in sorted(cnt) if k <= 6})
