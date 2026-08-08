#!/usr/bin/env python3
"""run_persist.py -- exact persistency on small modmuls + the conditioned curve.

Everything is gated: no fix is reported unless it holds in every enumerated
ground state."""
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import persist
from collections import defaultdict, Counter


KW = dict(mult='schoolbook', leaf=8, red='naf', mode='wallace')


def kindcount(Q, fix):
    c = Counter(Q.kind[v] for v in fix)
    return dict(c)


def one_prime(p):
    Q, A, B, C, xs = persist.ground_states(p, **KW)
    nstates = len(xs)
    # unconditional exact ceiling
    ceil, _ = persist.ceiling(Q, xs)
    persist.verify_subset(ceil, Q, xs, f"ceiling p={p}")
    # constructive
    prop = persist.prop_fix(Q)
    persist.verify_subset(prop, Q, xs, f"prop p={p}")
    prob = persist.probe_fix(Q)
    persist.verify_subset(prob, Q, xs, f"probe p={p}")
    return dict(p=p, n=Q.n, states=nstates,
                ceil=len(ceil), prop=len(prop), probe=len(prob),
                ceil_kinds=kindcount(Q, ceil),
                probe_kinds=kindcount(Q, prob),
                A=A, B=B, xs=xs, Q=Q, ceilmap=ceil)


def conditioned(res, kbits):
    """ceiling conditioned on pinning the low k bits of operand a (A) to their
    value in a canonical ground state -- averaged over the distinct patterns
    present in the solution set."""
    Q, A, xs = res['Q'], res['A'], res['xs']
    s = len(A.bits)
    out = []
    for k in kbits:
        if k > s:
            continue
        # group ground states by the low-k pattern of A, take mean fixed count
        counts = []
        patterns = set()
        for x in xs:
            pat = tuple(x[A.bits[t]] for t in range(k))
            patterns.add(pat)
        for pat in patterns:
            restrict = {A.bits[t]: pat[t] for t in range(k)}
            fix, m = persist.ceiling(Q, xs, restrict=restrict)
            persist.verify_subset(fix, Q,
                                  [x for x in xs if all(x[v] == vv for v, vv in restrict.items())],
                                  f"cond p={res['p']} k={k}")
            counts.append(len(fix))
        out.append((k, sum(counts) / len(counts), min(counts), max(counts),
                    len(patterns)))
    return out


def main():
    primes = [13, 29, 61, 127, 251]
    print("EXACT PERSISTENCY on small modmuls (schoolbook/naf/wallace)")
    print(f"{'p':>5} {'vars':>6} {'#gnd':>7} {'ceiling':>8} {'prop':>6} "
          f"{'probe':>6}   ceiling-kinds")
    results = []
    for p in primes:
        try:
            r = one_prime(p)
        except AssertionError as e:
            print(f"  p={p}: {e}")
            raise
        results.append(r)
        print(f"{p:>5} {r['n']:>6} {r['states']:>7} {r['ceil']:>8} "
              f"{r['prop']:>6} {r['probe']:>6}   {r['ceil_kinds']}")

    print("\nCONDITIONED CEILING vs # pinned operand-a bits (exact, per prime)")
    print("  (mean fixed vars over the distinct low-k patterns; k=s means a fully pinned)")
    cond_dump = {}
    for r in results:
        curve = conditioned(r, [0, 1, 2, 4, r['Q'] and len(r['A'].bits)])
        cond_dump[r['p']] = curve
        s = len(r['A'].bits)
        cells = "  ".join(f"k={k}:{mean:.0f}[{lo}-{hi}]" for k, mean, lo, hi, np in curve)
        print(f"  p={r['p']:>4} (s={s}, {r['n']} vars): {cells}")

    dump = [{kk: vv for kk, vv in r.items()
             if kk in ('p', 'n', 'states', 'ceil', 'prop', 'probe',
                       'ceil_kinds', 'probe_kinds')} for r in results]
    with open(os.path.join(os.path.dirname(__file__), 'persist_small.json'), 'w') as f:
        json.dump({'unconditional': dump,
                   'conditioned': {str(k): v for k, v in cond_dump.items()}}, f, indent=2)
    return results


if __name__ == '__main__':
    sys.setrecursionlimit(1000000)
    main()
