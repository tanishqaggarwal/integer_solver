#!/usr/bin/env python3
"""Add the HIDDEN linear constraints (Q's of the 72 (linear)^2 atoms) to the linear
system and re-solve over GF(P). Check: does the rank increase? Do the twist vars get
determined in terms of control bits? Does any constraint directly link the 22-side
(x_9770,x_3183) to the 233-side (x_18274,x_17728)?"""
import json, time
from collections import defaultdict
from propagate import load_atoms, atom_vars, NVARS
from check_square import try_sqrt

P = (1 << 61) - 1
def inv(a): return pow(a % P, P-2, P)

def main():
    t0 = time.time()
    atoms = load_atoms()
    rows = []
    # purely-linear atoms
    for poly in atoms:
        if all(len(m) <= 1 for m in poly):
            row = {}; const = 0
            for m, c in poly.items():
                if len(m) == 0: const = (const+c) % P
                else: row[m[0]] = (row.get(m[0], 0)+c) % P
            row = {v: c for v, c in row.items() if c % P}
            if row: rows.append((row, (-const) % P))
    npure = len(rows)
    # hidden linear constraints from (linear)^2 atoms
    nhidden = 0
    hidden_rows = []
    for poly in atoms:
        deg = max(len(m) for m in poly) if poly else 0
        if deg == 2 and len(poly) >= 3:
            Q = try_sqrt(poly)
            if Q and max(len(m) for m in Q) == 1:
                row = {}; const = 0
                for m, c in Q.items():
                    if len(m) == 0: const = (const+c) % P
                    else: row[m[0]] = (row.get(m[0], 0)+c) % P
                row = {v: c for v, c in row.items() if c % P}
                if row:
                    hidden_rows.append((row, (-const) % P)); nhidden += 1
    print(f"purely-linear rows {npure}, hidden-linear rows {nhidden} ({time.time()-t0:.0f}s)", flush=True)

    def gauss(all_rows):
        pivot_of = {}
        for row, rhs in all_rows:
            r = dict(row); rr = rhs; changed = True
            while changed:
                changed = False
                for v in list(r):
                    if v in pivot_of:
                        f = r.pop(v); prow, prhs = pivot_of[v]
                        for k, cv in prow.items():
                            nv = (r.get(k, 0) - f*cv) % P
                            if nv: r[k] = nv
                            elif k in r: del r[k]
                        rr = (rr - f*prhs) % P; changed = True; break
            if not r:
                continue
            pv = min(r); f = inv(r[pv])
            pivot_of[pv] = ({k: (cv*f) % P for k, cv in r.items() if k != pv}, (rr*f) % P)
        return pivot_of

    piv_before = gauss(rows)
    piv_after = gauss(rows + hidden_rows)
    print(f"rank before {len(piv_before)}, after adding hidden {len(piv_after)} (delta {len(piv_after)-len(piv_before)}) ({time.time()-t0:.0f}s)", flush=True)

    control = set(json.load(open('control_bits.json')))
    BITS22 = {1782,1858,2795,2800,3483,5443,10652,19520,21188,21588,23634,26947,27512,29682,30104,30596,30658,30792,33251,37748,37885,38116}
    # check if twist vars are determined and their control-bit support after
    import sys; sys.setrecursionlimit(200000)
    memo = {}
    def reduce_var(v):
        if v in memo: return memo[v]
        if v not in piv_after: memo[v] = ({v: 1}, 0); return memo[v]
        prow, prhs = piv_after[v]
        acc = defaultdict(int); const = prhs
        for k, c in prow.items():
            sub, sc = reduce_var(k)
            for fv, fc in sub.items(): acc[fv] = (acc[fv] + c*fc) % P
            const = (const + c*sc) % P
        memo[v] = ({k: c for k, c in acc.items() if c % P}, const % P); return memo[v]
    for w in (9770, 3183, 18274, 17728):
        f, c = reduce_var(w)
        cb = [v for v in f if v in control]
        c22 = [v for v in cb if v in BITS22]; c233 = [v for v in cb if v not in BITS22]
        # is it linked to the other side's vars?
        print(f"  x_{w}: {len(f)} free vars; control bits {len(cb)} ({len(c22)}x22, {len(c233)}x233)", flush=True)
    # KEY: does x_9770 - x_18274 reduce to something small (a twist-linking constraint)?
    f9, c9 = reduce_var(9770); f18, c18 = reduce_var(18274)
    gap = defaultdict(int)
    for k, v in f9.items(): gap[k] = (gap[k]+v) % P
    for k, v in f18.items(): gap[k] = (gap[k]-v) % P
    gap = {k: v for k, v in gap.items() if v % P}
    print(f"  x_9770 - x_18274 reduces to {len(gap)} free vars (const {(c9-c18)%P})", flush=True)
    f3, c3 = reduce_var(3183); f17, c17 = reduce_var(17728)
    gap2 = defaultdict(int)
    for k, v in f3.items(): gap2[k] = (gap2[k]+v) % P
    for k, v in f17.items(): gap2[k] = (gap2[k]-v) % P
    gap2 = {k: v for k, v in gap2.items() if v % P}
    print(f"  x_3183 - x_17728 reduces to {len(gap2)} free vars (const {(c3-c17)%P})", flush=True)
    print(f"done ({time.time()-t0:.0f}s)", flush=True)

if __name__ == '__main__':
    main()
