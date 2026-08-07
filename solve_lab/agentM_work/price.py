"""CANDIDATE-AGNOSTIC PRICING PRIMITIVE.

Input: a set of HANDLE VARIABLES to corrupt.  No assumption about where they come from
-- they need not touch the current failing equations.  Output: an exactly-measured score.

Why the caller only has to name handles
---------------------------------------
Corrupting handle h means demoting its definer atom so h becomes free.  Any OTHER variable
whose definer atom references h would otherwise absorb the corruption (its value would be
rewritten), so it must be demoted too and its value chosen rather than derived.  That
collateral set is computed here from the definer-level user graph -- the caller does not
supply it.  For the deliverable's four handles it adds exactly one variable, x_7068,
which is precisely engine2's fifth demotion, DERIVED rather than assumed.

Pricing paths
-------------
price_given(handles, values)  : values supplied -> forward, re-propagate, exact score.
price_search(handles)         : values NOT supplied -> seed the freed vars at their
                                uncorrupted values (all demoted atoms zero, so the state
                                is the uncorrupted baseline), then solve in equation space
                                over subsets of the live failures, APPLY, re-propagate and
                                measure.  Collateral damage is measured, never assumed.
"""
import sys, os, json, time, math, itertools, collections, pickle
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentE_work')
sys.set_int_max_str_digits(20_000_000)
import harness as H
import engine as EB
import engine3 as E3
import fscore, sparse

NEQ = len(H.eqt)
NV = H.NV

# definer-level users: w is a user of u if w's DEFINER atom mentions u
_USERS = collections.defaultdict(list)
for _w in H.SEQ:
    _i, _k = H.definer[_w]
    for _u in H.avars[_i]:
        if _u != _w:
            _USERS[_u].append(_w)

VAR2ATOM = {u: H.definer[u][0] for u in H.SEQ}


def closure(handles, depth=1, cap=64):
    """Handles + the variables that would absorb their corruption, to `depth` levels.

    depth=0 : handles only -- corruption propagates freely downstream.
    depth=1 : handles + their DIRECT definer-level users.  This is the containment the
              39,026 point uses: for its four handles the only direct user is x_7068,
              so depth=1 reproduces engine2's demotion set exactly.
    Deeper levels are available but over-demote: the 39,026 point lets its 18 further
    downstream variables propagate normally (their definer atoms stay satisfied), so
    demoting them is unnecessary, and for the value-search path it changes the baseline.
    """
    freed = set(handles)
    frontier = set(handles)
    for _ in range(max(0, depth)):
        nxt = set()
        for u in frontier:
            for w in _USERS.get(u, ()):
                if w not in freed and w in VAR2ATOM:
                    nxt.add(w)
        if not nxt:
            break
        freed |= nxt
        frontier = nxt
        if len(freed) > cap:
            return None, None            # corruption not containable within cap
    demote = sorted(VAR2ATOM[u] for u in freed if u in VAR2ATOM)
    return sorted(freed), demote


class Pricer:
    def __init__(self, base_vector):
        self.vd = base_vector
        # uncorrupted baseline: E's ORIGINAL orientation, every definer atom zero
        self.seed_unc = {f: base_vector[f] for f in EB.FREE if base_vector[f] != 0}
        self.v_unc = EB.forward(self.seed_unc)

    # ---------- path 1: values supplied ----------
    def price_given(self, handles, values, extra_seed=None, depth=1):
        freed, demote = closure(handles, depth=depth)
        if freed is None:
            return {'ok': False, 'why': 'closure exceeded cap'}
        eng = E3.Eng(demote)
        seed = {f: self.vd[f] for f in eng.FREE if self.vd[f] != 0}
        if extra_seed:
            seed.update(extra_seed)
        for u, val in values.items():
            if val:
                seed[u] = val
            else:
                seed.pop(u, None)
        v = eng.forward(seed)
        bad = eng.badatoms(v)
        fl = sorted(fscore.fails(bad))
        return {'ok': True, 'freed': freed, 'demote': demote, 'score': NEQ - len(fl),
                'fails': fl, 'nbad': len(bad), 'bad': sorted(bad), 'v': v, 'eng': eng,
                'seed': seed}

    # ---------- path 2: values searched ----------
    def price_search(self, handles, maxk=2, budget=90.0, want_vec=False, depth=1):
        t0 = time.time()
        freed, demote = closure(handles, depth=depth)
        if freed is None:
            return {'ok': False, 'why': 'closure exceeded cap'}
        eng = E3.Eng(demote)
        # seed freed vars at their UNCORRUPTED values -> all demoted atoms are zero
        seed = {f: self.vd[f] for f in eng.FREE if self.vd[f] != 0}
        for u in freed:
            if self.v_unc[u]:
                seed[u] = self.v_unc[u]
            else:
                seed.pop(u, None)
        v0 = eng.forward(seed)
        bad0 = eng.badatoms(v0)
        FAILS = sorted(fscore.fails(bad0))
        base_sc = NEQ - len(FAILS)
        if not FAILS:
            return {'ok': True, 'score': base_sc, 'freed': freed, 'note': 'no failures',
                    'base_score': base_sc, 'secs': time.time() - t0}

        cand = set(freed)
        for e in FAILS:
            for c, a in H.eqt[e][2]:
                if a >= 0:
                    try:
                        cand |= set(EB.cone(a)[1])
                    except Exception:
                        pass
        FS = set(eng.FREE)
        cand = sorted(f for f in cand if f in FS)
        aff, cols = _affine_cols(eng, v0, bad0, cand)

        CMx = {}
        for e in FAILS:
            cm = collections.defaultdict(int); const = 0
            for c, a in H.eqt[e][2]:
                if a < 0:
                    const += c
                else:
                    cm[a] += c
            CMx[e] = (dict(cm), const)

        def row_for(e):
            cm, const = CMx[e]
            row = {}
            for f in aff:
                co = 0
                for a, d in cols[f].items():
                    c = cm.get(a)
                    if c:
                        co += c * d
                if co:
                    row[f] = co
            s0 = const + sum(c * bad0[a] for a, c in cm.items() if a in bad0)
            return row, -s0

        ROWS = {e: row_for(e) for e in FAILS}
        best = (base_sc, None, None)
        for k in range(1, min(maxk, len(FAILS)) + 1):
            for S in itertools.combinations(FAILS, k):
                if time.time() - t0 > budget:
                    break
                rows = [ROWS[e][0] for e in S]; rhs = [ROWS[e][1] for e in S]
                if any(not r for r in rows):
                    continue
                sol, msg, _ = sparse.solve_sparse(rows, rhs, verbose=False,
                                                  maxcore=400, maxcorebits=5_000_000)
                if sol is None:
                    continue
                ns = dict(seed)
                for f, d in sol.items():
                    if d:
                        ns[f] = v0[f] + d
                try:
                    v = eng.forward(ns)
                    av = eng.badatoms(v)
                    sc = fscore.score(av)
                except Exception:
                    continue
                if sc > best[0]:
                    best = (sc, S, ns)
            if time.time() - t0 > budget:
                break
        out = {'ok': True, 'score': best[0], 'via': best[1], 'base_score': base_sc,
               'freed': freed, 'demote': demote, 'nknobs': len(aff),
               'secs': time.time() - t0}
        if want_vec and best[2] is not None:
            out['seed'] = best[2]
            out['eng'] = eng
        return out


def _resid_delta(eng, v0, base_bad, changes):
    v, aff = eng.apply_delta(v0, changes)
    touched = set()
    for u in set(aff) | set(changes):
        touched.update(E3._ATOM_OF[u])
    ns = {'v': v, '__builtins__': {}}
    bad = dict(base_bad)
    for i in touched:
        r = eval(H.acodes[i], ns)
        if r:
            bad[i] = r
        else:
            bad.pop(i, None)
    return bad, v


def _affine_cols(eng, v0, bad0, cand):
    cols = {}; aff = []
    for f in cand:
        o = v0[f]
        try:
            b1, _ = _resid_delta(eng, v0, bad0, {f: o + 1})
            b2, _ = _resid_delta(eng, v0, bad0, {f: o + 2})
            b7, _ = _resid_delta(eng, v0, bad0, {f: o + 7})
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
    return aff, cols


def load_deliverable():
    d = json.load(open('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json'))
    v = [0] * NV
    for k, val in d.items():
        v[int(k.split('_')[1])] = int(val)
    return v
