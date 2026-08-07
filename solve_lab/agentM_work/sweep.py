"""Exhaustive pricing of all C(32,4)=35,960 four-handle sites from the incident pool.

Speed work (validated against the slow path before use):
  * scoring a greedy prefix used a full engine3.forward + badatoms (~0.5 s).  Only the
    freed variables move, so an incremental resid_delta gives the same bad-atom dict for
    a fraction of the cost.  Verified identical on the calibration site.
  * engine3.Eng rebuilt its definer-level user map per site (30k-var loop).  The map does
    not depend on the demotion set except that pinned variables cannot be targets, so it
    is built ONCE here and pinned targets are skipped during traversal.

Ordering: 4-subsets containing at least one of the deliverable's four are priced FIRST, so
an interrupted run still yields a meaningful ordered prefix with a stated stopping point.
"""
import sys, os, json, time, math, itertools, collections, pickle
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentE_work')
sys.set_int_max_str_digits(20_000_000)
import harness as H
import engine as EB, engine3 as E3
import price as PR, fscore, sparse

NEQ = PR.NEQ
D4 = [642, 28730, 29854, 31864]

# ---- global definer-level user map, built once ----
BASE_USERS = collections.defaultdict(list)
for _w in H.SEQ:
    _i, _k = H.definer[_w]
    for _u in H.avars[_i]:
        if _u != _w:
            BASE_USERS[_u].append(_w)
BASE_USERS = dict(BASE_USERS)


def fast_downstream(eng, changed, pinset):
    aff = set(); st = list(changed)
    while st:
        u = st.pop()
        for w in BASE_USERS.get(u, ()):
            if w not in aff and w not in pinset:
                aff.add(w); st.append(w)
    return aff


def fast_apply(eng, v0, changes, pinset):
    v = list(v0)
    for k, val in changes.items():
        v[k] = val
    aff = fast_downstream(eng, changes.keys(), pinset)
    ns = {'v': v, '__builtins__': {}}
    for u in sorted(aff, key=lambda u: eng._pos[u]):
        i, kind = eng.definer[u]
        E3._solvevar(v, ns, u, i, kind[0])
    return v, aff


def fast_resid(eng, v0, bad0, changes, pinset):
    v, aff = fast_apply(eng, v0, changes, pinset)
    touched = set()
    for u in set(aff) | set(changes):
        touched.update(E3._ATOM_OF[u])
    ns = {'v': v, '__builtins__': {}}
    bad = dict(bad0)
    for i in touched:
        r = eval(H.acodes[i], ns)
        if r:
            bad[i] = r
        else:
            bad.pop(i, None)
    return bad, v


class Ctx:
    """Everything that does not depend on the site."""
    def __init__(self):
        self.vd = PR.load_deliverable()
        self.seed_unc = {f: self.vd[f] for f in EB.FREE if self.vd[f] != 0}
        self.v_unc = EB.forward(self.seed_unc)


def tune_fast(ctx, handles, nprobe=10, budget=60.0, want=False):
    t0 = time.time()
    freed, demote = PR.closure(handles)
    if freed is None:
        return {'ok': False, 'why': 'closure cap'}
    eng = E3.Eng(demote)
    pinset = set(eng.pin)
    seed = {f: ctx.vd[f] for f in eng.FREE if ctx.vd[f] != 0}
    for u in freed:
        if ctx.v_unc[u]:
            seed[u] = ctx.v_unc[u]
        else:
            seed.pop(u, None)
    v0 = eng.forward(seed)
    bad0 = eng.badatoms(v0)
    FAILS = sorted(fscore.fails(bad0))
    base = NEQ - len(FAILS)

    # affine columns of the freed handles only
    cols = {}; aff = []
    for f in sorted(freed):
        o = v0[f]
        try:
            b1, _ = fast_resid(eng, v0, bad0, {f: o + 1}, pinset)
            b2, _ = fast_resid(eng, v0, bad0, {f: o + 2}, pinset)
            b7, _ = fast_resid(eng, v0, bad0, {f: o + 7}, pinset)
        except Exception:
            continue
        col = {}; ok = True
        for a in set(b1) | set(b2) | set(b7) | set(bad0):
            d1 = b1.get(a, 0) - bad0.get(a, 0)
            if b2.get(a, 0) - bad0.get(a, 0) != 2 * d1 or b7.get(a, 0) - bad0.get(a, 0) != 7 * d1:
                ok = False; break
            if d1:
                col[a] = d1
        if ok:
            aff.append(f); cols[f] = col

    rows, rhs = [], []
    for e in FAILS:
        cm = collections.defaultdict(int); const = 0
        for c, a in H.eqt[e][2]:
            if a < 0:
                const += c
            else:
                cm[a] += c
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
        rhs.append(-(const + sum(c * bad0[a] for a, c in cm.items() if a in bad0)))

    order = [i for i in range(len(rows)) if rows[i]]
    if not order:
        return {'ok': True, 'score': base, 'base_score': base, 'nrows_target': 0,
                'nknobs': len(aff), 'secs': time.time() - t0, 'freed': freed}

    keep, sols = [], []
    for i in order:
        if time.time() - t0 > budget * 0.7:
            break
        trial = keep + [i]
        s, _, _ = sparse.solve_sparse([rows[j] for j in trial], [rhs[j] for j in trial],
                                      verbose=False, maxcore=400, maxcorebits=5_000_000)
        if s is not None:
            keep = trial
            sols.append(s)

    best = (base, None, None)
    if sols:
        idx = sorted(set([len(sols) - 1] +
                         [round(k * (len(sols) - 1) / max(1, nprobe - 1)) for k in range(nprobe)]))
        for j in idx:
            ch = {}
            for f, d in sols[j].items():
                if d:
                    ch[f] = v0[f] + d
            if not ch:
                continue
            try:
                bad, _ = fast_resid(eng, v0, bad0, ch, pinset)
                sc = fscore.score(bad)
            except Exception:
                continue
            if sc > best[0]:
                best = (sc, j, ch)
            if time.time() - t0 > budget:
                break
    out = {'ok': True, 'score': best[0], 'base_score': base, 'nrows_target': len(order),
           'nknobs': len(aff), 'nsol': len(sols), 'secs': time.time() - t0, 'freed': freed}
    if want and best[2] is not None:
        ns = dict(seed); ns.update(best[2])
        out['seed'] = ns; out['eng'] = eng
    return out


def ordered_subsets(pool):
    d4 = set(D4)
    first, rest = [], []
    for c in itertools.combinations(sorted(pool), 4):
        (first if d4 & set(c) else rest).append(c)
    return first, rest


if __name__ == '__main__':
    mode = sys.argv[1] if len(sys.argv) > 1 else 'validate'
    pool = json.load(open('incident_pool.json'))['incident_handles']
    ctx = Ctx()

    if mode == 'validate':
        print('=== validate fast tuner against the slow path on the calibration site ===')
        t0 = time.time(); r = tune_fast(ctx, D4, want=True); dt = time.time() - t0
        print(f'  fast tuner: base {r["base_score"]} -> {r["score"]}  ({dt:.2f}s, '
              f'{r["nrows_target"]} rows, {r["nknobs"]} knobs)')
        print(f'  CALIBRATION {"PASSED" if r["score"] >= 39026 else "FAILED"}')
        # confirm the fast score equals a full re-propagation
        if r.get('seed') is not None:
            eng = r['eng']; v = eng.forward(r['seed'])
            full = fscore.score(eng.badatoms(v))
            print(f'  full re-propagation of the same seed: {full}  '
                  f'-> incremental scoring exact: {full == r["score"]}')
        first, rest = ordered_subsets(pool)
        print(f'\n  ordering: {len(first)} subsets touch the deliverable\'s four, '
              f'{len(rest)} do not, total {len(first)+len(rest)}')
        print('\n=== timing sample (20 sites) ===')
        ts = []
        for c in (first[:10] + rest[:10]):
            t0 = time.time(); tune_fast(ctx, list(c)); ts.append(time.time() - t0)
        avg = sum(ts) / len(ts)
        tot = avg * (len(first) + len(rest))
        print(f'  mean {avg:.2f}s/site  -> {tot/3600:.1f} core-hours, '
              f'{tot/3600/4:.1f} h on 4 cores')
