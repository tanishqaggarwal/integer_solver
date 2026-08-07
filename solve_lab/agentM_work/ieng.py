"""Incremental site engine — no per-site Eng build, no per-site full forward.

What made the general enumeration slow: every 4-subset changes the demotion set, so
`E3.Eng(demote)` was rebuilt (30k-entry SEQ/SOLVE/_pos) and a full `forward` (30,378
_solvevar calls) plus a full `badatoms` (40,727 evals) ran per site.

Two observations remove all of that:

1. **The baseline vector is the SAME for every site.**  Seeding the freed variables at their
   uncorrupted values and propagating reproduces `v_unc` exactly, whatever the demotion set,
   because `v_unc` satisfies every definer atom.  So `v_unc` and `badatoms(v_unc)` are
   computed ONCE and shared.  It follows that the baseline failing set (the 25) is shared too,
   so the equation coefficient maps are precomputed once as well.

2. **No engine object is needed.**  A site is fully described by its pinned set; propagation
   is the global `H.SEQ` order with pinned variables skipped (they are inputs, not solved).

Per site the cost is then: a few incremental probes for the affine columns, the small greedy
solves, and a few incremental scorings.  Nothing scales with 30k or 40k.
"""
import sys, os, json, math, time, collections
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentE_work')
sys.set_int_max_str_digits(20_000_000)
import harness as H
import engine as EB
import price as PR, fscore, sparse

NV = H.NV
NEQ = PR.NEQ
acodes = H.acodes
definer = H.definer

BASE_POS = {u: i for i, u in enumerate(H.SEQ)}
BASE_USERS = collections.defaultdict(list)
for _w in H.SEQ:
    _i, _k = definer[_w]
    for _u in H.avars[_i]:
        if _u != _w:
            BASE_USERS[_u].append(_w)
BASE_USERS = dict(BASE_USERS)

ATOM_OF = collections.defaultdict(list)
for _i, _vs in enumerate(H.avars):
    for _u in _vs:
        ATOM_OF[_u].append(_i)


def _solvevar(v, ns, u, i, kd):
    c = acodes[i]
    v[u] = 0; c0 = eval(c, ns)
    v[u] = 1; c1 = eval(c, ns)
    if kd == 'lin':
        sl = c1 - c0
        v[u] = -c0 // sl if sl and c0 % sl == 0 else 0
    else:
        v[u] = 2; c2 = eval(c, ns)
        A2 = c2 - 2 * c1 + c0; A = A2 // 2; B = c1 - c0 - A; C = c0
        disc = B * B - 4 * A * C
        if disc < 0 or A == 0:
            v[u] = 0; return
        r = math.isqrt(disc)
        if r * r != disc:
            v[u] = 0; return
        rts = {(-B + s) // (2 * A) for s in (r, -r) if (-B + s) % (2 * A) == 0}
        v[u] = rts.pop() if len(rts) == 1 else 0


# ---------- shared baseline ----------
VD = PR.load_deliverable()
SEED_UNC = {f: VD[f] for f in EB.FREE if VD[f] != 0}
V_UNC = EB.forward(SEED_UNC)
BAD_UNC = {}
_ns = {'v': V_UNC, '__builtins__': {}}
for _i in range(len(H.atoms)):
    _r = eval(acodes[_i], _ns)
    if _r:
        BAD_UNC[_i] = _r
FAILS_UNC = sorted(fscore.fails(BAD_UNC))

# equation coefficient maps for the shared baseline failures
CM = {}
for _e in FAILS_UNC:
    _cm = collections.defaultdict(int); _c0 = 0
    for _c, _a in H.eqt[_e][2]:
        if _a < 0:
            _c0 += _c
        else:
            _cm[_a] += _c
    CM[_e] = (dict(_cm), _c0)


def apply_delta(v0, changes, pinset):
    v = list(v0)
    for k, val in changes.items():
        v[k] = val
    aff = set(); st = list(changes)
    while st:
        u = st.pop()
        for w in BASE_USERS.get(u, ()):
            if w not in aff and w not in pinset:
                aff.add(w); st.append(w)
    ns = {'v': v, '__builtins__': {}}
    for u in sorted(aff, key=lambda x: BASE_POS[x]):
        i, kind = definer[u]
        _solvevar(v, ns, u, i, kind[0])
    return v, aff


def resid(v0, bad0, changes, pinset):
    v, aff = apply_delta(v0, changes, pinset)
    touched = set()
    for u in set(aff) | set(changes):
        touched.update(ATOM_OF[u])
    ns = {'v': v, '__builtins__': {}}
    bad = dict(bad0)
    for i in touched:
        r = eval(acodes[i], ns)
        if r:
            bad[i] = r
        else:
            bad.pop(i, None)
    return bad, v


def site(handles):
    """freed set (handles + depth-1 definer-level users) and the pinned set."""
    freed, demote = PR.closure(handles)
    if freed is None:
        return None, None
    return freed, set(freed)


def score_from_unc(changes, pinset):
    bad, v = resid(V_UNC, BAD_UNC, changes, pinset)
    return fscore.score(bad), bad, v


def price_given(handles, values):
    """Values supplied (the calibration path)."""
    freed, pin = site(handles)
    if freed is None:
        return None
    ch = {u: values[u] for u in freed if u in values}
    sc, bad, v = score_from_unc(ch, pin)
    return {'score': sc, 'fails': sorted(fscore.fails(bad)), 'bad': sorted(bad),
            'v': v, 'freed': freed}


def affine_cols(pin, freed):
    cols = {}; aff = []
    for f in sorted(freed):
        o = V_UNC[f]
        try:
            b1, _ = resid(V_UNC, BAD_UNC, {f: o + 1}, pin)
            b2, _ = resid(V_UNC, BAD_UNC, {f: o + 2}, pin)
            b7, _ = resid(V_UNC, BAD_UNC, {f: o + 7}, pin)
        except Exception:
            continue
        col = {}; ok = True
        for a in set(b1) | set(b2) | set(b7) | set(BAD_UNC):
            d1 = b1.get(a, 0) - BAD_UNC.get(a, 0)
            if b2.get(a, 0) - BAD_UNC.get(a, 0) != 2 * d1 or \
               b7.get(a, 0) - BAD_UNC.get(a, 0) != 7 * d1:
                ok = False; break
            if d1:
                col[a] = d1
        if ok:
            aff.append(f); cols[f] = col
    return aff, cols


def tune(handles, nprobe=10, budget=30.0, want=False):
    t0 = time.time()
    freed, pin = site(handles)
    if freed is None:
        return {'ok': False, 'why': 'closure cap'}
    aff, cols = affine_cols(pin, freed)
    base = NEQ - len(FAILS_UNC)
    if not aff:
        return {'ok': True, 'score': base, 'base_score': base, 'nrows_target': 0,
                'nknobs': 0, 'freed': freed, 'secs': time.time() - t0}
    rows = []; rhs = []
    for e in FAILS_UNC:
        cm, const = CM[e]
        row = {}
        for f in aff:
            co = 0
            for a, d in cols[f].items():
                c = cm.get(a)
                if c:
                    co += c * d
            if co:
                row[f] = co
        rows.append(row)
        rhs.append(-(const + sum(c * BAD_UNC[a] for a, c in cm.items() if a in BAD_UNC)))
    order = [i for i in range(len(rows)) if rows[i]]
    if not order:
        return {'ok': True, 'score': base, 'base_score': base, 'nrows_target': 0,
                'nknobs': len(aff), 'freed': freed, 'secs': time.time() - t0}
    keep = []; sols = []
    for i in order:
        if time.time() - t0 > budget * 0.7:
            break
        trial = keep + [i]
        s, _, _ = sparse.solve_sparse([rows[j] for j in trial], [rhs[j] for j in trial],
                                      verbose=False, maxcore=400, maxcorebits=5_000_000)
        if s is not None:
            keep = trial; sols.append(s)
    best = (base, None, None)
    if sols:
        idx = sorted(set([len(sols) - 1] +
                         [round(k * (len(sols) - 1) / max(1, nprobe - 1))
                          for k in range(nprobe)]))
        for j in idx:
            ch = {f: V_UNC[f] + d for f, d in sols[j].items() if d}
            if not ch:
                continue
            try:
                bad, _ = resid(V_UNC, BAD_UNC, ch, pin)
                sc = fscore.score(bad)
            except Exception:
                continue
            if sc > best[0]:
                best = (sc, j, ch)
            if time.time() - t0 > budget:
                break
    out = {'ok': True, 'score': best[0], 'base_score': base, 'nrows_target': len(order),
           'nknobs': len(aff), 'freed': freed, 'secs': time.time() - t0}
    if want and best[2] is not None:
        out['changes'] = best[2]; out['pin'] = pin
    return out
