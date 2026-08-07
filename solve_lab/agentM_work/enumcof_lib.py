"""The cofactor-widened pricer, factored out of enumcof.py so other scripts can import it.

Identical arithmetic to enumcof.py's `tune_ext`: knobs = ieng.site(W) closure u the cofactors,
affinity of every knob tested by second differences (ieng.affine_cols), greedy row keeping over
exact integer solves, score by re-propagation.
"""
import sys, os, json, time

MDIR = '/home/user/integer_solver/solve_lab/agentM_work'
sys.path.insert(0, MDIR)
import shim                                                    # noqa: F401
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentE_work')
sys.set_int_max_str_digits(20_000_000)
import ieng, fscore, sparse                                    # noqa: E402

COF12 = [105, 1329, 3387, 5081, 5676, 9413, 10903, 11436, 14393, 14768, 17325, 22820]
COF4 = [1329, 9413, 10903, 17325]
EXTRA = COF4

V_UNC, BAD_UNC, CM, FAILS_UNC = ieng.V_UNC, ieng.BAD_UNC, ieng.CM, ieng.FAILS_UNC
BASE = ieng.NEQ - len(FAILS_UNC)

PF = json.load(open(os.path.join(MDIR, 'pfamily.json')))
SETS = {k: sorted({v['h'] for v in PF[f'incident_{k2}'].values()})
        for k, k2 in (('12', '7'), ('16', '12'), ('18', '25'))}
HL = SETS['16']
D4 = [642, 28730, 29854, 31864]


def tune_ext(handles, nprobe=80, budget=180.0, extra=None):
    t0 = time.time()
    ex = EXTRA if extra is None else extra
    freed, pin = ieng.site(handles)
    if freed is None:
        return None
    knobs = sorted(set(freed) | set(ex))
    aff, cols = ieng.affine_cols(pin, knobs)
    if not aff:
        return {'score': BASE, 'changes': None, 'pin': pin, 'naff': 0}
    rows, rhs = [], []
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
        return {'score': BASE, 'changes': None, 'pin': pin, 'naff': len(aff)}
    keep, sols = [], []
    for i in order:
        if time.time() - t0 > budget * 0.7:
            break
        trial = keep + [i]
        s, _, _ = sparse.solve_sparse([rows[j] for j in trial], [rhs[j] for j in trial],
                                      verbose=False, maxcore=400, maxcorebits=5_000_000)
        if s is not None:
            keep = trial; sols.append(s)
    best = (BASE, None)
    if sols:
        Lm = len(sols) - 1
        idx = sorted(set([Lm] + [round(k * Lm / max(1, nprobe - 1)) for k in range(nprobe)]))
        for j in idx:
            ch = {f: V_UNC[f] + d for f, d in sols[j].items() if d}
            if not ch:
                continue
            try:
                bad, _ = ieng.resid(V_UNC, BAD_UNC, ch, pin)
                sc = fscore.score(bad)
            except Exception:
                continue
            if sc > best[0]:
                best = (sc, ch)
            if time.time() - t0 > budget:
                break
    return {'score': best[0], 'changes': best[1], 'pin': pin, 'naff': len(aff)}
