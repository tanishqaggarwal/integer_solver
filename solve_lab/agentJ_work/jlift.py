#!/usr/bin/env python3
"""Integer lift: given free-input integer values, propagate over Z and solve the
free 'handle' variables so every constraint holds exactly in Z."""
import os, pickle, sys, time
from collections import defaultdict, deque
import jengine as E, jman as J, jmodp as MP

P = MP.P
NV = E.NV
definer, order, FREE = J.definer, J.order, J.FREE
polys = E.polys
CONS = MP.CONS

pos = {v: k for k, v in enumerate(order)}
uses = defaultdict(list)          # var -> defined vars that read it
for v, i in definer.items():
    for w in E.varsof[i]:
        if w != v:
            uses[w].append(v)
occ = defaultdict(list)
for i, s in enumerate(E.varsof):
    for v in s:
        occ[v].append(i)
EV = J.ev


def fwd_from(val, seeds):
    """recompute only the forward cone of `seeds` (topological)."""
    dirty = set()
    q = deque(seeds)
    while q:
        x = q.popleft()
        for w in uses[x]:
            if w not in dirty and w not in FREE:
                dirty.add(w)
                q.append(w)
    for v in sorted(dirty, key=lambda z: pos[z]):
        c, rest = EV[v]
        s = 0
        for k, cc in rest:
            t = cc
            for j in k:
                t *= val[j]
            s += t
        val[v] = -s // c if c in (1, -1) else (-s) // c
    return val


def atomval(i, val):
    s = 0
    for k, c in polys[i].items():
        t = c
        for j in k:
            t *= val[j]
        s += t
    return s


# --- structural handle map: constraint -> list of candidate free handle vars ---
def cone_free(i, maxdepth=5):
    seen = set(); out = []
    q = deque((w, 0) for w in E.varsof[i])
    while q:
        x, d = q.popleft()
        if x in seen: continue
        seen.add(x)
        if x in FREE:
            if len(occ[x]) == 1:
                out.append(x)
            continue
        if d >= maxdepth: continue
        di = definer.get(x)
        if di is None: continue
        for w in E.varsof[di]:
            if w != x and w not in seen:
                q.append((w, d + 1))
    return out


HCACHE = {}
def handles(i):
    if i not in HCACHE:
        HCACHE[i] = cone_free(i)
    return HCACHE[i]


def lift(val, rounds=6, verbose=True):
    val = list(val)
    E.forward(val, order, J.ev, definer, FREE)
    for r in range(rounds):
        bad = [i for i in CONS if atomval(i, val) != 0]
        if verbose:
            print(f"  lift round {r}: {len(bad)} constraints nonzero over Z")
        if not bad:
            break
        fixed = 0
        for i in bad:
            v0 = atomval(i, val)
            if v0 == 0:
                continue
            for z in handles(i):
                old = val[z]
                val[z] = old + 1
                fwd_from(val, [z])
                v1 = atomval(i, val)
                d = v1 - v0
                if d != 0 and v0 % d == 0:
                    val[z] = old - v0 // d
                    fwd_from(val, [z])
                    if atomval(i, val) == 0:
                        fixed += 1
                        break
                val[z] = old
                fwd_from(val, [z])
        if verbose:
            print(f"     fixed {fixed}")
        if fixed == 0:
            break
    return val


if __name__ == '__main__':
    import jsolve2
    b1, b2 = (int(sys.argv[1]), int(sys.argv[2])) if len(sys.argv) > 2 else (1, 0)
    mval, bad = jsolve2.branch(b1, b2)
    print("mod-p violated:", bad)
    val = list(J.BASE)
    for z in FREE:
        val[z] = mval[z] % P if mval[z] else 0
    t0 = time.time()
    val = lift(val)
    s, fails, av = E.score(val)
    print(f"score {s}/39033  fails={len(fails)}  ({time.time()-t0:.0f}s)")
    print("fails:", fails[:30])
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), f'J_b{b1}{b2}_{s}.json')
    E.save(val, out)
    print("saved", out)
