"""Corrected channel model + simultaneous solve, on engine2 (represents the deliverable)."""
import sys, os, json, re, collections, itertools, time, pickle
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentE_work')
sys.set_int_max_str_digits(20_000_000)
import harness as H
import engine as E_broken          # only for cone() (orientation-independent: uses definer)
import engine2 as E2, fast2, sparse
import mcore2 as M

P = M.P
ROWS = M.ROWS
CLUSTERKN = [6083, 11436, 14393, 14853, 22820, 26489, 31339, 37012]
# the 5 vars engine2 promotes to free are prime knobs: they drive the 8 residual atoms
PINKN = list(E2.PIN)


def channels(seed, atomset=None):
    """Channel classes at this seed (corrected engine). ON leaves reported separately."""
    v0, bad0, sig = M.measure(seed, atomset=atomset)
    cls = collections.defaultdict(list)
    on = []
    inert = []
    for f, c in sig.items():
        if c == 'ON':
            on.append(f); continue
        if any(c):
            cls[c].append(f)
        else:
            inert.append(f)
    return v0, bad0, cls, on, inert


def affine_cols(v0, bad0, cand):
    """Knobs whose effect on the residual is exactly affine, with their columns."""
    cols = {}; aff = []
    for f in cand:
        o = v0[f]
        b1, _ = fast2.resid_delta(v0, bad0, {f: o + 1})
        b2, _ = fast2.resid_delta(v0, bad0, {f: o + 2})
        b7, _ = fast2.resid_delta(v0, bad0, {f: o + 7})
        col = {}; ok = True
        for a in set(b1) | set(b2) | set(b7) | set(bad0):
            d1 = b1.get(a, 0) - bad0.get(a, 0)
            if b2.get(a, 0) - bad0.get(a, 0) != 2 * d1 or b7.get(a, 0) - bad0.get(a, 0) != 7 * d1:
                ok = False; break
            if d1:
                col[a] = d1
        if ok:
            aff.append(f); cols[f] = col
    return aff, cols


def simsolve(seed, maxr=3, maxv=2000, extrakn=()):
    """Simultaneous linear repair of the residual, on the corrected engine."""
    v0 = E2.forward(seed)
    bad0 = E2.badatoms(v0)
    if not bad0:
        return 0, seed, [], v0
    seedkn = list(CLUSTERKN) + list(PINKN) + list(extrakn)
    S = set(seedkn); pend = set(bad0); seenA = set(); cols = {}; knobs = []
    for rnd in range(maxr + 1):
        new = set()
        for a in pend:
            new |= set(E_broken.cone(a)[1])
        new -= S | set(seed)
        new = {f for f in new if f in E2.FREE or E2.definer[f] is None}
        if not new:
            break
        aff, c2 = affine_cols(v0, bad0, sorted(new))
        cols.update(c2); knobs += aff; S |= set(new)
        t = set()
        for f in aff:
            t |= set(cols[f])
        seenA |= pend
        pend = (t | set(bad0)) - seenA
        if len(S) > maxv:
            break
    aff0, c0 = affine_cols(v0, bad0, seedkn)
    for f in aff0:
        if f not in cols:
            knobs.append(f); cols[f] = c0[f]
    rows_at = set(bad0)
    for f in knobs:
        rows_at |= set(cols[f])
    rows_at = sorted(rows_at)
    rowmap = {a: {} for a in rows_at}
    for f in knobs:
        for a, c in cols[f].items():
            rowmap[a][f] = c
    sol, msg, _ = sparse.solve_sparse([rowmap[a] for a in rows_at],
                                      [-bad0.get(a, 0) for a in rows_at],
                                      names=rows_at, verbose=False,
                                      maxcore=400, maxcorebits=5_000_000)
    if sol is None:
        keep = []
        for i, a in enumerate(rows_at):
            idx = keep + [i]
            s2, _, _ = sparse.solve_sparse([rowmap[rows_at[j]] for j in idx],
                                           [-bad0.get(rows_at[j], 0) for j in idx],
                                           verbose=False, maxcore=400, maxcorebits=5_000_000)
            if s2 is not None:
                keep = idx
        if not keep:
            return None
        sol, _, _ = sparse.solve_sparse([rowmap[rows_at[j]] for j in keep],
                                        [-bad0.get(rows_at[j], 0) for j in keep],
                                        verbose=False, maxcore=400, maxcorebits=5_000_000)
    if sol is None:
        return None
    ns = dict(seed)
    for f, d in sol.items():
        if d:
            ns[f] = v0[f] + d
    v = E2.forward(ns)
    av = E2.badatoms(v)
    return len(E2.eqfails(av)), ns, sorted(av), v


def dump(v, path):
    json.dump({f"x_{i}": int(v[i]) for i in range(E2.NV) if v[i] != 0}, open(path, 'w'))
