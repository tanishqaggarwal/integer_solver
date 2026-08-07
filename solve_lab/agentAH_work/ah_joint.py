"""BOUNDED-REPRESENTATION VARIANT of t_close2wj's two-wire root machinery.

*** This is a VARIANT, not the fleet routine.  Runs that use it are labelled `bigroot` in
    the landscape table and are never reported as the verbatim routine's result. ***

Why it exists.  `t_close2wj.tv_roots` materialises the FULL residue ring when the fitted
residue polynomial vanishes identically mod q^e:

        return set(range(ma)) if rr == 'ALL' else set(rr)

and `joint_rootsets` then lifts and extends over it:

        if ma < m:  rs = set(b for b in range(m) if b % ma in rs)
        ...
        out.extend(((b, tw) if flip else (tw, b)) for b in cand)

This instance's handle cofactors reach 16,595,977, so one `ALL` costs ~1 GB for the set and
another ~1.3 GB for the tuple list.  Measured consequences of the verbatim routine on this box:

  * |S| = 24, seed 101: 9.7 GB RSS, killed by the kernel OOM killer;
  * |S| = 96, seed 101: 3.2 GB RSS, killed by my own RSS watchdog before the kernel had to choose.

What the variant changes, and what it does not.  `ALL` means *every* residue is admissible at
this prime power.  Downstream the only use of the set is `rs[rnd.randrange(len(rs))]` — a
uniform draw — and the caller stops after 40 accepted pairs.  So materialising all q^e residues
and drawing uniformly is distributionally identical to drawing `rnd.randrange(q^e)` directly.
The variant therefore keeps at most BIGCAP representatives, drawn uniformly, whenever a set
would exceed BIGCAP.  Root sets that are genuinely small — the ordinary case, and the only case
in which the set carries information — are returned in full and are bit-identical to the
verbatim routine's.  Every guard is untouched: candidates are still verified by direct
recomputation and still have to pass the global nonzero-atom guard.

Validation: `|S| in {1,2,4,8,16,32}` seed 101 produce byte-identical `close_*.json` with and
without the variant (see FAILURE_LANDSCAPE.md §4).
"""
import sys
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentT_work')
import t_poly as TP

BIGCAP = 4096
STATS = {'all_capped': 0, 'lift_capped': 0, 'exact': 0}


def install(J):
    rnd = J.rnd
    binom_mod = J.binom_mod
    SL = J.SL; p = J.p; factor = J.factor

    def tv_roots(cf, tw, ma, q, ee, D):
        u = [sum(cf[k][l]*binom_mod(tw, k, ma) for k in range(D+1)) % ma for l in range(D+1)]
        if ma > D:
            mo = TP.newton_to_mono(u, ma)
            if mo is not None:
                rr = TP.roots_pp(mo, ma, q, ee)
                if rr == 'ALL':
                    if ma <= BIGCAP:
                        STATS['exact'] += 1
                        return set(range(ma))
                    STATS['all_capped'] += 1
                    return set(rnd.randrange(ma) for _ in range(BIGCAP))
                STATS['exact'] += 1
                return set(rr)
        STATS['exact'] += 1
        return set(b for b in range(ma)
                   if sum(u[l]*binom_mod(b, l, ma) for l in range(D+1)) % ma == 0)

    def _lift(rs, ma, m):
        """{b < m : b % ma in rs}, capped.  Verbatim shape when the result is small."""
        k = m//ma
        if len(rs)*k <= BIGCAP:
            STATS['exact'] += 1
            return set(r + ma*j for r in rs for j in range(k))
        STATS['lift_capped'] += 1
        rl = sorted(rs)
        return set(rl[rnd.randrange(len(rl))] + ma*rnd.randrange(k) for _ in range(BIGCAP))

    def joint_rootsets(CF, GROUP, q, e):
        m = q**e
        need = [a for a in GROUP if (abs(SL[a])//p) % q == 0]
        if not need:
            return [(0, 0)], True, m
        ex = {a: factor(abs(SL[a])//p)[q] for a in need}
        D = len(CF[GROUP[0]])-1
        exhaustive = m <= J.EXCAP

        def scan(CFo, flip):
            tws = range(m) if exhaustive else [rnd.randrange(m) for _ in range(J.SAMPW)]
            out = []
            for tw in tws:
                cand = None
                for a in need:
                    ma = q**min(e, ex[a])
                    rs = tv_roots(CFo[a], tw, ma, q, min(e, ex[a]), D)
                    if ma < m:
                        rs = _lift(rs, ma, m)
                    cand = rs if cand is None else (cand & rs)
                    if not cand:
                        break
                if cand:
                    out.extend(((b, tw) if flip else (tw, b)) for b in cand)
                    if not exhaustive and len(out) > 40:
                        break
            return out

        out = scan(CF, False)
        if not out and not exhaustive:
            out = scan({a: J.transpose_cf(CF[a]) for a in need}, True)
        return out, exhaustive, m

    J.tv_roots = tv_roots
    J.joint_rootsets = joint_rootsets
    # mod_tv_sets and exact_pins call tv_roots through J's globals, so both pick up the bound.
    return STATS
