"""jm step 9: build an EXPLICIT sparse move direction over F_p and measure it.

Iterative column growth (frame2.py's loop, but solving for the KERNEL):
  columns C = support(target);   rows R = all checks except the allowed-to-break
  set D.   Solve  A_{R,C} c = 0,  T_C c = 1.   If infeasible, add the columns
  occurring in the rows that block it, and repeat.
Then apply the direction exactly, fwd2, repair, and count equations.

CHUNKED/RESUMABLE: every attempt is appended to jm_constr.jsonl.
usage: python3 jm_09_construct.py <target C1|C2> [maxcols]
"""
import os, sys, json, time, collections, random
sys.path.insert(0, '/home/user/integer_solver/solve_lab/s10')
import jm_lib as J
import lib as L, tools as T_, ad
from jm_07_lin import load, weight, FREEATOMS
import jm_05_engine as EN
P = J.P
LOG = '/home/user/integer_solver/solve_lab/s10/jm_constr.jsonl'


def gauss(rows, ncols):
    """rows: list of (dict col->coef, rhs).  Returns (solution dict, None) or
    (None, blocking row indices)."""
    piv = {}
    order = []
    blocked = []
    for idx, (r, rhs) in enumerate(rows):
        r = {k: v % P for k, v in r.items() if v % P}
        rhs = rhs % P
        while True:
            hit = None
            for c in r:
                if c in piv:
                    hit = c
                    break
            if hit is None:
                break
            f = r[hit]
            pr, prh = piv[hit]
            for k, v in pr.items():
                nv = (r.get(k, 0) - f * v) % P
                if nv:
                    r[k] = nv
                else:
                    r.pop(k, None)
            rhs = (rhs - f * prh) % P
        if not r:
            if rhs:
                blocked.append(idx)
            continue
        c0 = min(r)
        inv = pow(r[c0], -1, P)
        piv[c0] = ({k: v * inv % P for k, v in r.items()}, rhs * inv % P)
        order.append(c0)
    if blocked:
        return None, blocked
    sol = {}
    for c0 in reversed(order):
        pr, prh = piv[c0]
        s = prh
        for k, v in pr.items():
            if k != c0:
                s -= v * sol.get(k, 0)
        sol[c0] = s % P
    return {k: v for k, v in sol.items() if v % P}, None


def build_direction(target, D, rows, cols, U, maxit=25, seed_cols=None,
                    verbose=True):
    """returns (c dict u->val, predicted broken rows) or (None, blockers)"""
    Tv = target
    C = set(seed_cols if seed_cols is not None else Tv)
    R = [a for a in rows if a not in D]
    for it in range(maxit):
        Cl = sorted(C)
        sysrows = []
        rowid = []
        for a in R:
            rr = {u: v for u, v in rows[a].items() if u in C}
            if rr:
                sysrows.append((rr, 0))
                rowid.append(a)
        sysrows.append(({u: v for u, v in Tv.items() if u in C}, 1))
        rowid.append('T')
        sol, blk = gauss(sysrows, len(Cl))
        if sol is not None:
            return sol, None
        bad = [rowid[i] for i in blk]
        if verbose:
            print(f'   it{it}: {len(C)} cols, blocked by {bad[:6]}'
                  f'{"..." if len(bad) > 6 else ""}', flush=True)
        # grow C by the neighbourhood: every row touching C contributes its cols
        new = set()
        for a in bad:
            if a != 'T':
                new |= set(rows[a])
        for a in R:
            if new and len(new - C) > 300:
                break
            if C & set(rows[a]):
                new |= set(rows[a])
        new -= C
        if not new:
            return None, bad
        C |= new
    return None, ['maxit']


def measure(c, W0, tag, extra=None):
    v = list(W0)
    for u, d in c.items():
        v[u] = v[u] + d
    J.fwd2(v, 2)
    c1, s1, nz1, av1 = EN.state(v)
    rec = {'tag': tag, 'support': len(c), 'raw_out12': c1, 'raw_score': s1,
           'raw_broken': nz1}
    if extra:
        rec.update(extra)
    print(f'  [{tag}] support={len(c)} raw out12={c1} score={s1} broken={nz1[:12]}',
          flush=True)
    return v, rec


if __name__ == '__main__':
    which = sys.argv[1] if len(sys.argv) > 1 else 'C1'
    t0 = time.time()
    cols, dR, dA1 = load()
    U = sorted(cols)
    rows = collections.defaultdict(dict)
    for u, cc in cols.items():
        for a, d in cc.items():
            rows[a][u] = d
    rows = {a: r for a, r in rows.items() if a not in FREEATOMS}
    Wt = {a: weight(a) for a in rows}
    Tv = ({u: dR[u][0] for u in U if dR[u][0]} if which == 'C1'
          else {u: dA1[u] for u in U if dA1[u]})
    print(f'target {which}: {sorted(Tv)}', flush=True)
    W0 = J.base_state()
    f = open(LOG, 'a')

    Dsets = [([41400, 41507, 41827, 42245], 'greedy4'),
             ([41400, 41507, 41827, 42245, 3576], 'greedy4+3576'),
             ([41400, 41507, 41827, 42245, 7930, 41512], 'greedy4+7930'),
             ([a for a in rows if Wt[a] <= 1], 'all-weight1')]
    for D, nm in Dsets:
        print(f'\n--- D = {nm} ({len(D)} rows, weight '
              f'{sum(Wt[a] for a in D if a in Wt)}) ---', flush=True)
        c, blk = build_direction(Tv, set(D), rows, cols, U)
        if c is None:
            print(f'   INFEASIBLE (blocked by {blk[:6]})', flush=True)
            f.write(json.dumps({'target': which, 'D': nm,
                                'feasible': False, 'blk': [str(x) for x in blk[:12]]}) + '\n')
            f.flush()
            continue
        pred = [a for a in rows
                if sum(rows[a].get(u, 0) * v for u, v in c.items()) % P]
        pe = set()
        for a in pred:
            pe |= set(L.atom2eq.get(a, {})) - J.E12
        print(f'   FEASIBLE: support {len(c)}; predicted broken {pred} '
              f'-> {len(pe)} eqs', flush=True)
        v, rec = measure(c, W0, f'{which}-{nm}',
                         {'pred_rows': pred, 'pred_eqs': len(pe),
                          'c': {str(k): str(x) for k, x in c.items()}})
        keep = EN.keep_c1 if which == 'C1' else EN.keep_c2
        rec['keeps'] = bool(keep(v))
        vr, c2, s2, nz2 = EN.repair(list(v), keep=keep, tag=nm, verbose=False)
        rec['rep_out12'], rec['rep_score'], rec['rep_broken'] = c2, s2, nz2
        print(f'   repaired: out12={c2} score={s2} broken={nz2[:12]} '
              f'keeps={rec["keeps"]}  ({time.time()-t0:.0f}s)', flush=True)
        f.write(json.dumps(rec) + '\n')
        f.flush()
    f.close()
