"""Monotone closure: only accept a repair that strictly reduces the bad-check count."""
import sys, os, json, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L, fw, deep
P = L.P
NAT = {u: len(L.var_atoms[u]) for u in range(L.NVARS)}


def shallow(a, locked):
    out = []
    for u in L.avars[a]:
        if L.definer.get(u) is not None or u in locked:
            continue
        if any(mm.count(u) > 1 for mm in L.polys[a]):
            continue
        out.append(u)
    out.sort(key=lambda u: (NAT[u], u))
    return out


def close(v, locked, rounds=60, verbose=True, use_deep=True):
    t0 = time.time()
    cur = set(fw.bad_checks(v))
    for rnd in range(rounds):
        if not cur:
            break
        prog = False
        for a in sorted(cur, key=lambda a: (len(L.atom2eq.get(a, {})), a)):
            if fw.evalpoly(L.polys[a], v) == 0:
                continue
            cands = [(t, None) for t in shallow(a, locked)]
            if use_deep:
                try:
                    hs, base = deep.handles(v, a, locked=locked)
                    cands += [(t, d) for t, d in sorted(hs, key=lambda kv: (NAT[kv[0]], kv[0]))]
                except Exception:
                    pass
            for t, d in cands:
                old = v[t]
                if d is None:
                    x = fw.solve_lin(a, t, v)
                    if x is None or x == old:
                        continue
                else:
                    bs = fw.evalpoly(L.polys[a], v)
                    if d == 0 or bs % d:
                        continue
                    x = old - bs // d
                v[t] = x
                fw.forward(v)
                new = set(fw.bad_checks(v))
                if len(new) < len(cur):
                    cur = new
                    prog = True
                    break
                v[t] = old
                fw.forward(v)
        if verbose:
            print(f"  round{rnd}: bad={len(cur)} ({time.time()-t0:.0f}s) {sorted(cur)[:14]}", flush=True)
        if not prog:
            break
    fw.forward(v)
    return v, fw.bad_checks(v)
