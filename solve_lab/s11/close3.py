"""Staged closure with globally-unique handle assignment."""
import sys, os, json, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L, fw, deep
P = L.P
NAT = {u: len(L.var_atoms[u]) for u in range(L.NVARS)}


def cands(v, a, locked, used, use_deep=True):
    """(var, delta_or_None) candidates for zeroing atom a, cheapest first"""
    out = []
    for u in L.avars[a]:
        if L.definer.get(u) is not None or u in locked or u in used:
            continue
        if any(mm.count(u) > 1 for mm in L.polys[a]):
            continue
        out.append((u, None))
    out.sort(key=lambda kv: (NAT[kv[0]], kv[0]))
    if use_deep:
        try:
            hs, base = deep.handles(v, a, locked=set(locked) | set(used))
            out += [(t, d) for t, d in sorted(hs, key=lambda kv: (NAT[kv[0]], kv[0]))]
        except Exception:
            pass
    return out


def apply_cand(v, a, t, d):
    old = v[t]
    if d is None:
        x = fw.solve_lin(a, t, v)
        if x is None or x == old:
            return False
    else:
        bs = fw.evalpoly(L.polys[a], v)
        if d == 0 or bs % d:
            return False
        x = old - bs // d
    v[t] = x
    fw.forward(v)
    if fw.evalpoly(L.polys[a], v) == 0:
        return True
    v[t] = old
    fw.forward(v)
    return False


def close(v, locked, rounds=25, verbose=True):
    t0 = time.time()
    best = None
    for rnd in range(rounds):
        bad = fw.bad_checks(v)
        if not bad:
            break
        if best is None or len(bad) < best[0]:
            best = (len(bad), [x for x in v])
        used = set()
        # order: fewest candidate handles first (most constrained), then fewest equations
        order = sorted(bad, key=lambda a: (len(cands(v, a, locked, used, use_deep=False)),
                                           len(L.atom2eq.get(a, {})), a))
        for a in order:
            if fw.evalpoly(L.polys[a], v) == 0:
                continue
            for t, d in cands(v, a, locked, used):
                if apply_cand(v, a, t, d):
                    used.add(t)
                    break
        bad2 = fw.bad_checks(v)
        if verbose:
            print(f"  round{rnd}: bad={len(bad2)} ({time.time()-t0:.0f}s) {bad2[:14]}", flush=True)
        if len(bad2) == len(bad) and set(bad2) == set(bad):
            break
    fw.forward(v)
    b = fw.bad_checks(v)
    if best is not None and len(b) > best[0]:
        v[:] = best[1]
        fw.forward(v)
        b = fw.bad_checks(v)
    return v, b
