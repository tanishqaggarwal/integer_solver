"""Generic dependency-ordered repair driver with 1-step lookahead (movevar style)."""
import sys, collections, time
sys.path.insert(0, '/home/user/integer_solver/solve_lab/s9/quad')
from common import *

CODES = None


def score(v, tag='', quiet=False):
    nz = nz_checks(v); ng = nz_gates(v)
    live = [a for a in range(len(polys)) if evalpoly(polys[a], v) != 0]
    ff = H.evaluate(CODES, v, eqs_of(live))
    if not quiet:
        print(f'[{tag:34s}] nzcheck={len(nz):3d} nzgate={len(ng)} FAIL={len(ff):4d}  {sorted(set(nz+ng))[:24]}')
    return nz, ng, ff


FREEZE = set()


def absorbers(a, v):
    """free inputs appearing linearly in atom a's residual that can zero it."""
    R = resid_poly.get(a, polys[a])
    out = []
    for u in sorted(set(x for m in R for x in m)):
        if u not in freeset or u in FREEZE:
            continue
        c = 0; nl = False
        for m, cc in R.items():
            if len(m) == 1 and m[0] == u: c += cc
            elif u in m: nl = True
        if nl or c == 0: continue
        old = v[u]; v[u] = 0; rest = evalpoly(R, v); v[u] = old
        if rest % c: continue
        out.append((u, -rest // c))
    return out


def cone_absorbers(a, v, maxfree=400):
    """free inputs in the backward cone of atom a whose unit move changes it linearly."""
    R = resid_poly.get(a, polys[a])
    rts = sorted(set(x for m in R for x in m))
    seen = set(rts); q = collections.deque(rts); free = []
    while q:
        u = q.popleft()
        d = definer.get(u)
        if d is None:
            free.append(u); continue
        for w in set(x for m in polys[d] for x in m):
            if w != u and w not in seen:
                seen.add(w); q.append(w)
    r0 = evalpoly(R, v)
    out = []
    for u in free[:maxfree]:
        if u in FREEZE: continue
        vv = list(v); ripple(vv, {u: v[u] + 1})
        d1 = evalpoly(R, vv) - r0
        if d1 == 0: continue
        if r0 % d1 == 0:
            out.append((u, v[u] - r0 // d1, d1))
    return out


def repair(v, rounds=25, verbose=True, use_cone=True):
    hist = []
    for r in range(rounds):
        nz, ng, ff = score(v, f'round {r}', quiet=not verbose)
        hist.append((len(nz), len(ff)))
        if verbose:
            print(f'  round {r}: nz={len(nz)} FAIL={len(ff)}  {sorted(nz)[:20]}')
        if not nz: return v, hist
        best = None
        for a in nz:
            cands = [(u, val, None) for u, val in absorbers(a, v)]
            if use_cone and not cands:
                cands = cone_absorbers(a, v)
            for u, val, _ in cands:
                if v[u] == val: continue
                vv = list(v); ripple(vv, {u: val})
                n2 = len(nz_checks(vv)) + len(nz_gates(vv))
                if best is None or n2 < best[0]:
                    best = (n2, u, val, a)
        if best is None:
            if verbose: print('  no absorber found -> stall')
            return v, hist
        n2, u, val, a = best
        if n2 >= len(nz) + len(ng):
            if verbose: print(f'  best move x_{u} (atom {a}) gives {n2} >= {len(nz)+len(ng)} -> stall')
            return v, hist
        if verbose: print(f'  apply x_{u} := (for atom {a}) -> nz {n2}')
        ripple(v, {u: val})
    return v, hist
