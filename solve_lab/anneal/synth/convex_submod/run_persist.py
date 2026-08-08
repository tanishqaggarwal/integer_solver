#!/usr/bin/env python3
"""run_persist.py -- exact persistency on small modmuls + the conditioned curve.

Gated: no fix is reported unless it holds in every enumerated ground state."""
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import persist
from collections import Counter

KW = dict(mult='schoolbook', leaf=8, red='naf', mode='wallace')


def kindcount(Q, fix):
    return dict(Counter(Q.kind[v] for v in fix))


def crosscheck_p13():
    import verify
    Q, A, B, C = persist.make(13, **KW)
    s = (13).bit_length()
    fix, cnt = persist.ceiling(Q, A, B, C, 13, s)
    xs = verify.zero_states(Q)
    seen0 = bytearray(Q.n); seen1 = bytearray(Q.n)
    for x in xs:
        for v in range(Q.n):
            (seen1 if x[v] else seen0)[v] = 1
    dfs_fix = {v: (1 if seen1[v] else 0) for v in range(Q.n)
               if not (seen0[v] and seen1[v])}
    ok = (fix == dfs_fix)
    print(f"CROSS-CHECK p=13: replay-ceiling={len(fix)} dfs-ceiling={len(dfs_fix)} "
          f"replay-states={cnt} dfs-states={len(xs)} -> {'MATCH' if ok else 'MISMATCH'}")
    assert ok


def main():
    sys.setrecursionlimit(1000000)
    crosscheck_p13()
    primes = [13, 29, 61, 127, 251]
    print("\nEXACT PERSISTENCY on small modmuls (schoolbook/naf/wallace)")
    print(f"{'p':>5} {'s':>2} {'vars':>6} {'#gnd':>7} {'ceiling':>8} {'prop':>6} "
          f"{'probe':>6}   ceiling-kinds")
    results = []
    for p in primes:
        s = p.bit_length()
        Q, A, B, C = persist.make(p, **KW)
        fix, cnt = persist.ceiling(Q, A, B, C, p, s)
        persist.verify_subset(fix, Q, A, B, C, p, s, label=f"ceiling p={p}")
        prop = persist.prop_fix(Q)
        persist.verify_subset(prop, Q, A, B, C, p, s, label=f"prop p={p}")
        if Q.n <= 400:
            prob = persist.probe_fix(Q)
            persist.verify_subset(prob, Q, A, B, C, p, s, label=f"probe p={p}")
            probn = len(prob)
        else:
            probn = None
        results.append(dict(p=p, s=s, n=Q.n, states=cnt, ceil=len(fix),
                            prop=len(prop), probe=probn,
                            ceil_kinds=kindcount(Q, fix),
                            Q=Q, A=A, B=B, C=C))
        pv = probn if probn is not None else '-'
        print(f"{p:>5} {s:>2} {Q.n:>6} {cnt:>7} {len(fix):>8} {len(prop):>6} "
              f"{str(pv):>6}   {kindcount(Q, fix)}")

    print("\nCONDITIONED CEILING vs # pinned low bits of operand a (exact)")
    print("  fixed vars, averaged over the 2^k pin patterns; k=s => a fully known")
    cond_dump = {}
    for r in results:
        Q, A, B, C, s, p = r['Q'], r['A'], r['B'], r['C'], r['s'], r['p']
        ks = sorted(set([0, 1, 2, 4, s]))
        curve = []
        for k in ks:
            counts = []
            for pat in range(1 << k):
                pin = {t: (pat >> t) & 1 for t in range(k)}
                # exact by full enumeration; ceiling is the true constant-set
                fix, _ = persist.ceiling(Q, A, B, C, p, s, pin_a=pin)
                counts.append(len(fix))
            curve.append((k, sum(counts) / len(counts), min(counts), max(counts)))
        cond_dump[p] = curve
        cells = "  ".join(f"k={k}:{m:.0f}[{lo}-{hi}]" for k, m, lo, hi in curve)
        print(f"  p={p:>4} (s={s}, {Q.n} vars): {cells}")

    dump = [{kk: r[kk] for kk in ('p', 's', 'n', 'states', 'ceil', 'prop',
                                  'probe', 'ceil_kinds')} for r in results]
    with open(os.path.join(os.path.dirname(__file__), 'persist_small.json'), 'w') as f:
        json.dump({'unconditional': dump,
                   'conditioned': {str(k): v for k, v in cond_dump.items()}},
                  f, indent=2)


if __name__ == '__main__':
    main()
