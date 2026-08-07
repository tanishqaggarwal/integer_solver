"""Rank of the region response mod p, across STRUCTURALLY varied SELECTOR settings.

The one input held fixed in all 16 detach states is the selector setting: every configuration
priced so far inherits the deliverable's {2081, 24601}.  `p` enters the frame through its
constants, so this varies the setting structurally (position in the OR-tree recovered by
seltree.py, cardinality of the live set, which subtree it is drawn from) and re-measures the
region response.

Per configuration:

  region R      = every equation touching a nonzero atom of the state;
  b_i           = the exact value of region row i.  A row that is a single top-level SQUARE atom
                  is replaced by that atom's BASE -- `row = 0` iff `base = 0`.  Truncating the
                  square instead is the defect T caught in optN.inner and that I reintroduced in
                  pgap.price; it is handled here, not assumed away;
  M             = the EXACT linear part of those rows in the complete knob set (every free input
                  of Frame(POOL) syntactically supporting an atom of R).  Affineness is NOT
                  assumed: every (row, knob) entry is probed at t = 1, 2 and any entry whose
                  second difference is nonzero has its whole column re-read at t = 1..6 and the
                  degree-1 monomial coefficient recovered by EXACT Newton interpolation
                  (linear coeff = sum_k (-1)^(k-1) D^k g(0) / k), with integrality asserted.
                  A plain secant would be wrong in exactly those entries -- the same 28-of-7,399
                  error I found in widen.py.

Two levels:
  AMBIENT  -- all knobs free.  (At the deliverable this is already NOT rank-deficient mod p:
              rk_Q = rk_p = 12.  The deficiency is a property of the lattice, not the rows.)
  LATTICE  -- knobs restricted to the integer kernel of the LINEAR collateral response, i.e. the
              directions that provably keep every currently-satisfied outside equation satisfied
              to first order.  This is a RELAXATION of pgap.py's exact saturation lattice (it
              omits collateral constraints that are nonlinear at the base point), so it can only
              make the region look MORE solvable -- which is the safe direction for a hunt:
              gap_p = 1 on the relaxation implies gap_p >= 1 on the exact lattice.

Usage:  python3 pselrank.py <out.jsonl> [tag ...]      (skips tags already in out.jsonl)
"""
import os, sys, json, time
from collections import defaultdict
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
sys.set_int_max_str_digits(2000000)
import ev, optN
from optN import fr, FREE, FR0, BASEFV, atom_eqs, _bits
from frameB import State
import frameB as FB
from sqaudit import square_base
from ikc import int_kernel_columns
from flint import fmpz_mat
import re as _re
import psel

Pp = 115792089237316195423570985008687907853269984665640564039457584007908834671663
# CONTROL prime, same size, no relation to the instance: separates "the region response is
# deficient MOD p" from "the region response is simply low rank mod anything".
Qc = 115792089237316195423570985008687907853269984666640564039457584007908834671783
eq_terms = ev.eq_terms
_VR = _re.compile(r'x_(\d+)')
_BASECODE = {}
_MISS = object()
MAXT = 6            # probe points t = 1..MAXT; detects degree <= MAXT in one knob


def base_code(a):
    c = _BASECODE.get(a, _MISS)
    if c is _MISS:
        sb = square_base(a)
        _BASECODE[a] = c = (compile(_VR.sub(r'v[\1]', sb), '<sb>', 'eval') if sb else None)
    return c


def rank_mod(rows, ncol, p):
    A = [[x % p for x in r] for r in rows]
    r = 0
    for c in range(ncol):
        piv = None
        for i in range(r, len(A)):
            if A[i][c]:
                piv = i
                break
        if piv is None:
            continue
        A[r], A[piv] = A[piv], A[r]
        inv = pow(A[r][c], p - 2, p)
        A[r] = [(x * inv) % p for x in A[r]]
        for i in range(len(A)):
            if i != r and A[i][c]:
                f = A[i][c]
                A[i] = [(A[i][j] - f * A[r][j]) % p for j in range(ncol)]
        r += 1
        if r == len(A):
            break
    return r


def rank_q(rows, ncol):
    if not rows or not ncol:
        return 0
    return int(fmpz_mat([[int(x) for x in r] for r in rows]).rank())


def lin_coeff(vals):
    """vals = [g(0), g(1), ..., g(m)] of an integer polynomial g of degree <= m.
    Returns (c1, deg_ok) with c1 the EXACT monomial coefficient of t^1."""
    m = len(vals) - 1
    d = list(vals)
    diffs = [d[0]]
    for k in range(1, m + 1):
        d = [d[i + 1] - d[i] for i in range(len(d) - 1)]
        diffs.append(d[0])
    c1 = Fraction(0)
    for k in range(1, m + 1):
        if diffs[k]:
            c1 += Fraction((-1) ** (k - 1), k) * diffs[k]
    assert c1.denominator == 1, 'non-integer linear coefficient %s' % c1
    return int(c1), diffs[m] == 0


class Probe:
    """One-knob probe that skips State's equation bookkeeping: only the knob's DAG descendants
    and the atoms it can reach are recomputed."""

    def __init__(s, st, specR, specO):
        s.st = st
        s.specR = specR
        s.specO = specO
        s.baseav = st.av

    def read(s, Y, t):
        st = s.st
        v = st.v[:]
        ns = {'v': v, '__builtins__': {}}
        v[Y] = st.fv.get(Y, 0) + t
        for u in fr.desc[Y]:
            v[u] = eval(FB.DEFEXPR[u], ns)
        av = {a: eval(FB.ACODE[a], ns) for a in fr.chk.get(Y, [])}
        base = s.baseav
        out = []
        for specs in (s.specR, s.specO):
            row = []
            for rs in specs:
                if rs[1] is not None:
                    row.append(eval(rs[1], ns))
                else:
                    tot = 0
                    for c, a in eq_terms[rs[0]][2]:
                        x = av.get(a)
                        if x is None:
                            x = base.get(a)
                        if x:
                            tot += c * x
                    row.append(tot)
            out.append(row)
        return out


def rowspecs(st, Rl):
    """(equation, compiled base code or None) per row, plus a count of rooted rows."""
    out = []
    nroot = 0
    for e in Rl:
        m, sq, tl = eq_terms[e]
        live = [(c, a) for c, a in tl if st.av.get(a)]
        code = None
        if len(live) == 1 and base_code(live[0][1]) is not None:
            code = base_code(live[0][1])
            nroot += 1
        out.append((e, code))
    return out, nroot


def readrows0(st, specs):
    v = []
    for e, code in specs:
        if code is not None:
            v.append(eval(code, st.ns))
        else:
            t = 0
            for c, a in eq_terms[e][2]:
                x = st.av.get(a)
                if x:
                    t += c * x
            v.append(t)
    return v


MODE = ['pinned']       # 'pinned' = psel.state_for ; 'consistent' = pselc.state_for


def build_state(on):
    if MODE[0] == 'consistent':
        import pselc
        return pselc.state_for(on)
    return psel.state_for(on)


def measure(tag, on, lattice=True, cap_knobs=2400, cap_rows=3000, verbose=True):
    t0 = time.time()
    st = build_state(on)
    NZ, R = psel.region_of(st)
    Rl = sorted(R)
    if not Rl:
        return dict(tag=tag, note='empty region')
    if len(Rl) > cap_rows:
        return dict(tag=tag, nlive=len(on), R=len(Rl), note='region above cap_rows')
    K = psel.knobs_of(R)
    if len(K) > cap_knobs:
        return dict(tag=tag, nlive=len(on), R=len(Rl), knobs=len(K), note='knobs above cap')

    touched = set()
    for Y in K:
        for a in fr.chk.get(Y, []):
            touched |= atom_eqs[a]
    outside = sorted(touched - R)
    out_fail = len([e for e in st.fails if e not in R])

    specR, nroot = rowspecs(st, Rl)
    specO_all, _ = rowspecs(st, outside)
    b = readrows0(st, specR)
    bo_all = readrows0(st, specO_all)
    # only outside rows SATISFIED at the base state constrain the lattice
    keep = [i for i, x in enumerate(bo_all) if x == 0]
    specO = [specO_all[i] for i in keep]
    nR, nO = len(Rl), len(specO)
    bo = [0] * nO

    pr = Probe(st, specR, specO)
    cols, cols_o = [], []
    nonaff_entries = 0
    nonaff_rows = set()
    reinterp = 0
    deg_overflow = []
    for Y in K:
        r1, o1 = pr.read(Y, 1)
        r2, o2 = pr.read(Y, 2)
        cR = [r1[i] - b[i] for i in range(nR)]
        cO = [o1[i] - bo[i] for i in range(nO)]
        badR = [i for i in range(nR) if (r2[i] - b[i]) != 2 * cR[i]]
        badO = [i for i in range(nO) if (o2[i] - bo[i]) != 2 * cO[i]]
        if badR or badO:
            reinterp += 1
            nonaff_entries += len(badR) + len(badO)
            for i in badR:
                nonaff_rows.add(Rl[i])
            reads = [([b[i] for i in range(nR)], [bo[i] for i in range(nO)]), (r1, o1), (r2, o2)]
            for t in range(3, MAXT + 1):
                reads.append(pr.read(Y, t))
            for i in badR:
                c1, ok = lin_coeff([reads[t][0][i] for t in range(MAXT + 1)])
                cR[i] = c1
                if not ok:
                    deg_overflow.append((Y, Rl[i]))
            for i in badO:
                c1, ok = lin_coeff([reads[t][1][i] for t in range(MAXT + 1)])
                cO[i] = c1
                if not ok:
                    deg_overflow.append((Y, 'out%d' % i))
        cols.append(cR)
        cols_o.append(cO)

    n = len(K)
    M = [[cols[j][i] for j in range(n)] for i in range(nR)]
    rQ = rank_q(M, n)
    aug = [M[i] + [b[i]] for i in range(nR)]
    rQa = rank_q(aug, n + 1)
    rP = rank_mod(M, n, Pp)
    rPa = rank_mod(aug, n + 1, Pp)

    rec = dict(tag=tag, nlive=len(on), live=sorted(on), score=st.score(),
               R=nR, nz=len(NZ), knobs=n, rooted=nroot, outside=len(outside),
               collateral_live=nO, out_fail=out_fail, ceiling=39033 - out_fail,
               nonaffine_entries=nonaff_entries, nonaffine_knobs=reinterp,
               nonaffine_rows=sorted(nonaff_rows), degree_overflow=len(deg_overflow),
               amb_rk_Q=rQ, amb_rk_Q_aug=rQa, amb_gap_Q=rQa - rQ,
               amb_rk_p=rP, amb_rk_p_aug=rPa, amb_gap_p=rPa - rP,
               amb_deficiency=rQ - rP)

    if lattice:
        C = [[cols_o[j][i] for j in range(n)] for i in range(nO)]
        C = [r for r in C if any(r)]
        t1 = time.time()
        if C:
            Kn = int_kernel_columns(C, n)
        else:
            Kn = [[1 if i == j else 0 for i in range(n)] for j in range(n)]
        rec['lat_secs'] = round(time.time() - t1, 1)
        if not Kn:
            rec['lat_dim'] = 0
        else:
            Kb = fmpz_mat([[int(x) for x in v] for v in Kn]).lll().tolist()
            Kb = [[int(x) for x in r] for r in Kb if any(int(x) for x in r)]
            d = len(Kb)
            ML = [[sum(M[i][j] * Kb[a][j] for j in range(n)) for a in range(d)]
                  for i in range(nR)]
            augL = [ML[i] + [b[i]] for i in range(nR)]
            rec.update(lat_dim=d,
                       lat_rk_Q=rank_q(ML, d), lat_rk_Q_aug=rank_q(augL, d + 1),
                       lat_rk_p=rank_mod(ML, d, Pp), lat_rk_p_aug=rank_mod(augL, d + 1, Pp),
                       lat_rk_q=rank_mod(ML, d, Qc), lat_rk_q_aug=rank_mod(augL, d + 1, Qc))
            rec['lat_gap_q_ctl'] = rec['lat_rk_q_aug'] - rec['lat_rk_q']
            rec['lat_deficiency_ctl'] = rec['lat_rk_Q'] - rec['lat_rk_q']
            rec['lat_gap_Q'] = rec['lat_rk_Q_aug'] - rec['lat_rk_Q']
            rec['lat_gap_p'] = rec['lat_rk_p_aug'] - rec['lat_rk_p']
            rec['lat_deficiency'] = rec['lat_rk_Q'] - rec['lat_rk_p']
            rec['sat_rows'] = sum(1 for x in b if x == 0)
            # A region row whose ENTIRE lattice response vanishes but whose target does not can
            # never be zeroed on this lattice.  failing >= that count, so
            #   score <= 39033 - unzeroable.
            # Over Q this is exact; mod p it is the stronger test (M_i = 0 mod p, b_i != 0 mod p)
            # and is the row-level form of the p-obstruction.
            uq = sum(1 for i in range(nR) if b[i] != 0 and not any(ML[i]))
            up = sum(1 for i in range(nR)
                     if b[i] % Pp != 0 and not any(x % Pp for x in ML[i]))
            rec['unzeroable_Q'] = uq
            rec['unzeroable_p'] = up
            rec['score_ub_Q'] = 39033 - uq
            rec['score_ub_p'] = 39033 - up
    rec['secs'] = round(time.time() - t0, 1)
    if verbose:
        print(fmt(rec), flush=True)
    return rec


def fmt(r):
    if 'note' in r:
        return '%-24s SKIP %s' % (r['tag'], r['note'])
    s = ('%-24s live=%-3d score=%-6d |R|=%-4d k=%-4d na=%-4d | AMB rkQ=%-4d rkp=%-4d gapQ=%d '
         'gap_p=%d' % (r['tag'], r['nlive'], r['score'], r['R'], r['knobs'],
                       r['nonaffine_entries'], r['amb_rk_Q'], r['amb_rk_p'],
                       r['amb_gap_Q'], r['amb_gap_p']))
    if 'lat_dim' in r:
        s += ' | LAT dim=%-4d rkQ=%-4d rkp=%-4d gapQ=%s gap_p=%s ubp=%s' % (
            r['lat_dim'], r.get('lat_rk_Q', -1), r.get('lat_rk_p', -1),
            r.get('lat_gap_Q'), r.get('lat_gap_p'), r.get('score_ub_p'))
        s += ' | CTL rkq=%-4d gapq=%s' % (r.get('lat_rk_q', -1), r.get('lat_gap_q_ctl'))
    return s + ' (%.0fs)' % r.get('secs', 0)


def main():
    outp = os.path.join(HERE, 'runs', sys.argv[1] if len(sys.argv) > 1 else 'pselrank.jsonl')
    args = sys.argv[2:]
    if args and args[0] == 'consistent':
        MODE[0] = 'consistent'
        args = args[1:]
    shard = nshard = None
    if len(args) == 2 and args[0].isdigit() and args[1].isdigit():
        shard, nshard = int(args[0]), int(args[1])
        args = []
    want = set(args)
    done = set()
    # resumable: every completed tag already recorded IN THIS MODE is skipped
    pref = 'pselrankC' if MODE[0] == 'consistent' else 'pselrank_'
    for f in sorted(os.listdir(os.path.join(HERE, 'runs'))):
        if f.startswith(pref) and f.endswith('.jsonl'):
            for ln in open(os.path.join(HERE, 'runs', f)):
                try:
                    done.add(json.loads(ln)['tag'])
                except Exception:
                    pass
    fh = open(outp, 'a')
    # cheapest first, so a shard that is cut short has still answered the small configurations
    cfg = psel.configs()
    sizes = {r['tag']: r['knobs'] for r in json.load(open(os.path.join(HERE, 'runs',
                                                                      'psel_size.json')))}
    cfg.sort(key=lambda kv: sizes.get(kv[0], 10 ** 6))
    for i, (tag, on) in enumerate(cfg):
        if tag in done or (want and tag not in want):
            continue
        if nshard is not None and i % nshard != shard:
            continue
        try:
            rec = measure(tag, on)
        except Exception as ex:
            import traceback
            traceback.print_exc()
            rec = dict(tag=tag, note='ERROR %s: %s' % (type(ex).__name__, ex))
            print(fmt(rec), flush=True)
        fh.write(json.dumps(rec) + '\n')
        fh.flush()
    fh.close()


if __name__ == '__main__':
    main()
