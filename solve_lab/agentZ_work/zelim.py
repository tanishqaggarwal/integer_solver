#!/usr/bin/env python3
"""Agent Z, stage B: sparse Markowitz elimination of every NON-selector variable
from the booleanity-reduced linear system.  Whatever survives is a linear
constraint on the 256 selectors -- exactly where a cardinality constraint lives.

Strategy: repeatedly pick the non-selector column with the fewest incident rows,
and within it the sparsest row, as pivot.  Singleton columns (a variable defined
by exactly one row) are eliminated with zero fill, which is the common case here.
"""
import os, sys, json, pickle, collections, heapq

HERE = os.path.dirname(os.path.abspath(__file__))
sel = set(json.load(open(os.path.join(HERE, 'zsel.json')))['selectors'])
PRIME = (1 << 61) - 1
CONST = -1   # sentinel column for the constant term

def main():
    D = pickle.load(open(os.path.join(HERE, 'zlin.pkl'), 'rb'))
    lin_rows = D['lin_rows']
    print("linear rows:", len(lin_rows))

    rows = {}          # rid -> {col: coeff mod PRIME}
    colrows = collections.defaultdict(set)
    for rid, (eqi, q) in enumerate(lin_rows):
        r = {}
        for m, c in q.items():
            j = CONST if not m else m[0]
            r[j] = (r.get(j, 0) + c) % PRIME
        r = {j: c for j, c in r.items() if c}
        if not r:
            continue
        rows[rid] = r
        for j in r:
            if j != CONST:
                colrows[j].add(rid)
    print("rows:", len(rows), " columns:", len(colrows))
    nonsel_cols = set(colrows) - sel
    print("non-selector columns to eliminate:", len(nonsel_cols))

    heap = [(len(colrows[j]), j) for j in nonsel_cols]
    heapq.heapify(heap)
    elim = 0
    maxrowlen = 0
    while heap:
        n, j = heapq.heappop(heap)
        if j not in nonsel_cols:
            continue
        cur = colrows.get(j)
        if cur is None or len(cur) == 0:
            nonsel_cols.discard(j)
            continue
        if len(cur) != n:
            heapq.heappush(heap, (len(cur), j))
            continue
        # pivot row = sparsest row containing j
        prid = min(cur, key=lambda r: len(rows[r]))
        pr = rows[prid]
        inv = pow(pr[j], PRIME - 2, PRIME)
        pr = {k: (v * inv) % PRIME for k, v in pr.items()}
        rows[prid] = pr
        for rid in list(cur):
            if rid == prid:
                continue
            r = rows[rid]
            f = r.get(j, 0)
            if not f:
                colrows[j].discard(rid)
                continue
            for k, v in pr.items():
                old = r.get(k, 0)
                nv = (old - f * v) % PRIME
                if nv:
                    if not old:
                        r[k] = nv
                        if k != CONST:
                            colrows[k].add(rid)
                    else:
                        r[k] = nv
                else:
                    if old:
                        del r[k]
                        if k != CONST:
                            colrows[k].discard(rid)
            maxrowlen = max(maxrowlen, len(r))
        # retire the pivot row and column
        for k in pr:
            if k != CONST:
                colrows[k].discard(prid)
        del rows[prid]
        colrows.pop(j, None)
        nonsel_cols.discard(j)
        elim += 1
        if elim % 1000 == 0:
            print("   eliminated %d cols, rows left %d, maxrowlen %d"
                  % (elim, len(rows), maxrowlen), flush=True)

    print()
    print("=" * 78)
    print("non-selector columns eliminated:", elim)
    print("ROWS SURVIVING (support only in selectors + constant):", len(rows))
    nz = 0
    pure_const = 0
    empty0 = 0
    out = []
    for rid, r in rows.items():
        cols = sorted(k for k in r if k != CONST)
        assert all(c in sel for c in cols), ("leftover non-selector col", cols[:5])
        if not cols:
            pure_const += 1
            if not r:
                empty0 += 1
            continue
        nz += 1
        out.append((lin_rows[rid][0], sorted(r.items())))
    print("  survivors that are pure-constant rows:", pure_const, " of which EMPTY(0=0):", empty0, " of which INCONSISTENT(c=0,c!=0):", pure_const-empty0)
    print("  SURVIVORS THAT ARE GENUINE LINEAR CONSTRAINTS ON SELECTORS:", nz)
    for eqi, r in out[:60]:
        desc = []
        for k, v in r:
            vv = v if v < PRIME // 2 else v - PRIME
            desc.append(("%+d" % vv) if k == CONST else ("%+d*s%d" % (vv, k)))
        print("    eq%-6d : %s" % (eqi, " ".join(desc)[:400]))
    pickle.dump(out, open(os.path.join(HERE, 'zelim_survivors.pkl'), 'wb'))
    print("saved zelim_survivors.pkl")

if __name__ == '__main__':
    main()
