#!/usr/bin/env python3
"""Agent Z, TASK 3+4 decisive: solve the boolean layer SYMBOLICALLY.

The subsystem of equations whose variables are all boolean is (after x^2 -> x)
entirely LINEAR: 4390 rows over 2585 boolean variables, 256 of which are the
selectors.  Solve it treating the selectors as symbols, so every determined
boolean wire comes out as an affine function  a0 + sum a_i s_i.

Then:
  * a liveness/cardinality constraint would appear as a wire whose affine form
    has LARGE selector support (a count), or as a row with selector-only support.
  * a wire that is boolean-constrained but whose affine form can leave {0,1}
    would bound |S|.  We check every determined wire's affine form for that.
"""
import os, sys, json, pickle, collections, heapq
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
sel = set(json.load(open(os.path.join(HERE, 'zsel.json')))['selectors'])
SELL = sorted(sel)
SIDX = {s: i for i, s in enumerate(SELL)}

def main():
    D = pickle.load(open(os.path.join(HERE, 'zlin.pkl'), 'rb'))
    rows = [(i, q) for i, q in D['lin_rows'] if q]
    print("FULL booleanity-reduced linear rows:", len(rows))

    # matrix: unknown columns = non-selector boolean vars; RHS = affine in selectors
    M = {}      # rid -> ({col: Fraction}, rhs dict {selidx or -1 : Fraction})
    colrows = collections.defaultdict(set)
    for rid, (eqi, q) in enumerate(rows):
        lhs = {}
        rhs = {}
        for m, c in q.items():
            if not m:
                rhs[-1] = rhs.get(-1, Fraction(0)) - c          # move const to RHS
            elif m[0] in sel:
                rhs[SIDX[m[0]]] = rhs.get(SIDX[m[0]], Fraction(0)) - c
            else:
                lhs[m[0]] = lhs.get(m[0], Fraction(0)) + c
        lhs = {k: v for k, v in lhs.items() if v}
        rhs = {k: v for k, v in rhs.items() if v}
        M[rid] = (lhs, rhs)
        for k in lhs:
            colrows[k].add(rid)
    print("unknown wire columns:", len(colrows))

    heap = [(len(colrows[j]), j) for j in colrows]
    heapq.heapify(heap)
    pivot_of = {}
    live = set(colrows)
    elim = 0
    maxsup = 0
    while heap:
        n, j = heapq.heappop(heap)
        if j not in live:
            continue
        cur = colrows.get(j) or set()
        if not cur:
            live.discard(j); continue
        if len(cur) != n:
            heapq.heappush(heap, (len(cur), j)); continue
        prid = min(cur, key=lambda r: len(M[r][0]))
        plhs, prhs = M[prid]
        inv = Fraction(1) / plhs[j]
        plhs = {k: v * inv for k, v in plhs.items()}
        prhs = {k: v * inv for k, v in prhs.items()}
        M[prid] = (plhs, prhs)
        for rid in list(cur):
            if rid == prid:
                continue
            lhs, rhs = M[rid]
            f = lhs.get(j, 0)
            if not f:
                colrows[j].discard(rid); continue
            for k, v in plhs.items():
                old = lhs.get(k, Fraction(0))
                nv = old - f * v
                if nv:
                    lhs[k] = nv
                    if not old:
                        colrows[k].add(rid)
                elif old:
                    del lhs[k]; colrows[k].discard(rid)
            for k, v in prhs.items():
                nv = rhs.get(k, Fraction(0)) - f * v
                if nv:
                    rhs[k] = nv
                elif k in rhs:
                    del rhs[k]
            maxsup = max(maxsup, len(rhs))
        for k in plhs:
            colrows[k].discard(prid)
        pivot_of[j] = (prid, plhs, prhs)
        del M[prid]
        colrows.pop(j, None)
        live.discard(j)
        elim += 1
        if elim % 500 == 0:
            print("   pivots %d, rows left %d, max rhs support %d"
                  % (elim, len(M), maxsup), flush=True)

    print()
    print("pivots (wires solved):", elim, " rows left:", len(M))
    bad = 0
    selonly = 0
    for rid, (lhs, rhs) in M.items():
        assert not lhs, "leftover unknown"
        if rhs:
            selonly += 1
            print("  *** SELECTOR-ONLY SURVIVING ROW from eq", rows[rid][0], ":", sorted(rhs.items())[:20])
    print("SURVIVING ROWS WITH NONTRIVIAL SELECTOR SUPPORT:", selonly)

    # back-substitute pivots to fully affine forms
    order = list(pivot_of.keys())
    # iterate to fixpoint (pivot rows may reference later pivots)
    forms = {}
    changed = True
    rounds = 0
    pend = dict(pivot_of)
    while pend and rounds < 60:
        rounds += 1
        prog = 0
        for j in list(pend):
            prid, plhs, prhs = pend[j]
            others = [k for k in plhs if k != j]
            if all(k in forms for k in others):
                f = dict(prhs)
                for k in others:
                    c = plhs[k]
                    for t, v in forms[k].items():
                        nv = f.get(t, Fraction(0)) - c * v
                        if nv:
                            f[t] = nv
                        elif t in f:
                            del f[t]
                forms[j] = f
                del pend[j]
                prog += 1
        if not prog:
            break
    print("wires expressed as affine functions of the selectors:", len(forms),
          " (unresolved after", rounds, "rounds:", len(pend), ")")

    sup = collections.Counter(len([t for t in f if t != -1]) for f in forms.values())
    print("selector-support size of each solved wire:", sorted(sup.items()))
    # any wire that is a large weighted sum of selectors == a count
    big = [(j, f) for j, f in forms.items() if len([t for t in f if t != -1]) >= 3]
    print("wires depending affinely on >=3 selectors:", len(big))
    for j, f in big[:20]:
        print("   x%d = %s" % (j, " ".join(
            ("%+s" % v) if t == -1 else ("%+s*s%d" % (v, SELL[t])) for t, v in sorted(f.items()))[:200]))

    # boolean-range test: every solved wire must be in {0,1} for all sigma.
    viol = []
    for j, f in forms.items():
        const = f.get(-1, Fraction(0))
        coefs = [v for t, v in f.items() if t != -1]
        lo = const + sum(v for v in coefs if v < 0)
        hi = const + sum(v for v in coefs if v > 0)
        if lo < 0 or hi > 1 or any(v.denominator != 1 for v in coefs + [const]):
            viol.append((j, lo, hi, f))
    print()
    print("solved wires whose affine form can leave {0,1} over sigma in {0,1}^256:", len(viol))
    for j, lo, hi, f in viol[:20]:
        print("   x%d range [%s,%s] : %s" % (j, lo, hi, sorted(f.items())[:8]))
    pickle.dump({'forms': forms, 'pend': list(pend)}, open(os.path.join(HERE, 'zwire_full.pkl'), 'wb'))
    print("saved zwire.pkl")

if __name__ == '__main__':
    main()
