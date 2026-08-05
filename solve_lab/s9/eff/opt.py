"""Evaluate candidate sacrifice sets S: confined atoms -> knobs -> lattice -> max
simultaneously-zeroable subset (exact, over Z), with minimal-core pruning."""
import pickle, collections, itertools, sys, time
import lib as L, model as MD
from fractions import Fraction

v0 = None


def init():
    global v0
    if v0 is None:
        v0 = L.load(L.BEST24)
        MD.BASEP = [MD.prim_val(a, v0) for a in range(L.NA)]
    return v0


def rank_of(mod):
    A = mod['A']
    knobs = mod['knobs']
    if not knobs or not A:
        return 0
    m = [[Fraction(mod['M'][a].get(x, 0)) for x in knobs] for a in A]
    r = 0
    nr, nc = len(m), len(knobs)
    for c in range(nc):
        piv = next((i for i in range(r, nr) if m[i][c]), None)
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


def shrink(mod, T):
    T = set(T)
    for e in sorted(T):
        T2 = T - {e}
        if MD.solvable(mod, T2) is None:
            T = T2
    return frozenset(T)


def maximise(mod, S, floor=0, budget=400000, verbose=False):
    """Largest-first with minimal-unsolvable-core pruning."""
    S = sorted(S)
    n = len(S)
    bad = []
    calls = 0
    for size in range(n, floor - 1, -1):
        for T in itertools.combinations(S, size):
            Ts = set(T)
            if any(b <= Ts for b in bad):
                continue
            calls += 1
            if calls > budget:
                return None, None, None, calls          # gave up
            z = MD.solvable(mod, Ts)
            if z is not None:
                return size, T, z, calls
            core = shrink(mod, Ts)
            calls += len(Ts)
            bad.append(core)
            bad.sort(key=len)
            if len(bad) > 4000:
                bad = bad[:4000]
    return 0, (), {}, calls


FAIL0 = frozenset([9123, 9421, 12231, 12270, 12350, 14584, 18673, 22044, 29125])


def evaluate(S, label='', floor=0, budget=200000, verbose=True):
    init()
    S = frozenset(S)
    outside = len(FAIL0 - S)          # baseline-failing equations we did NOT sacrifice
    mod = MD.build(S, v0, verbose=False)
    D = rank_of(mod)
    sz, T, z, calls = maximise(mod, S, floor=floor, budget=budget)
    if sz is None:
        if verbose:
            print(f'{label} |S|={len(S)} |A|={len(mod["A"])} knobs={len(mod["knobs"])} D={D}  -> BUDGET EXCEEDED')
        return None
    fail = len(S) - sz + outside
    if verbose:
        print(f'{label} |S|={len(S)} |A|={len(mod["A"])} knobs={len(mod["knobs"])} D={D} '
              f'maxzero={sz} outside={outside} -> TOTAL FAILING={fail}  (calls={calls})')
    return {'S': S, 'mod': mod, 'D': D, 'maxzero': sz, 'T': T, 'z': z, 'fail': fail}
