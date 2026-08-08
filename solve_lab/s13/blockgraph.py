#!/usr/bin/env python3
"""
Answers two structural questions about the emitted QUBO set:

  Q1  "If the coupling is 80 bits, is the solve 2^80 operations?"
      The 80 was a COEFFICIENT MAGNITUDE (max |Q_ij| = 80, i.e. 7 bits of
      dynamic range), not a coupling width.  Coefficient size has NOTHING to do
      with search cost -- it is the analog precision an annealer must resolve.
      The quantity that does control exact solve cost is the SEPARATOR /
      TREEWIDTH of the block-interaction graph: exact dynamic programming over a
      tree decomposition costs ~2^(bits in the largest separator).
      This script MEASURES that, via min-degree and min-fill elimination
      (standard upper-bound heuristics for treewidth).

  Q2  "What is the size distribution of the QUBOs?"
      Full histogram of block sizes per tier.

Usage: python3 blockgraph.py
"""
import os, sys, json, time
from collections import defaultdict
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 's9', 'eff'))
import lib as L

NA = len(L.avars)
W = 16                       # bits per wire residue
DRIVERS = [2081, 4287, 8731, 9118]
D2 = [8731, 9118]


def downstream(seeds):
    users = defaultdict(list)
    for a in range(NA):
        for x in L.avars[a]:
            users[x].append(a)
    unknown, touched = set(seeds), set()
    dq = list(seeds)
    while dq:
        x = dq.pop()
        for a in users[x]:
            touched.add(a)
            oc = L.atom_out.get(a)
            if oc is not None and oc[1] not in unknown:
                unknown.add(oc[1]); dq.append(oc[1])
    return unknown, touched


def build_graph(atoms, unknown):
    """Nodes = atoms (blocks).  Edge = they share an UNKNOWN wire."""
    holders = defaultdict(set)
    for a in atoms:
        for x in L.avars[a]:
            if x in unknown:
                holders[x].add(a)
    adj = defaultdict(set)
    for x, ats in holders.items():
        ats = list(ats)
        for i in range(len(ats)):
            for j in range(i + 1, len(ats)):
                adj[ats[i]].add(ats[j]); adj[ats[j]].add(ats[i])
    for a in atoms:
        adj.setdefault(a, set())
    return adj, holders


def treewidth_ub(adj, method='min-fill'):
    """Upper bound on treewidth via elimination ordering (returns width in NODES)."""
    g = {v: set(ns) for v, ns in adj.items()}
    width = 0
    order = []
    while g:
        if method == 'min-degree':
            v = min(g, key=lambda u: len(g[u]))
        else:
            def fill(u):
                ns = list(g[u]); f = 0
                for i in range(len(ns)):
                    for j in range(i + 1, len(ns)):
                        if ns[j] not in g[ns[i]]:
                            f += 1
                return f
            v = min(g, key=fill)
        ns = g[v]
        width = max(width, len(ns))
        for a in ns:
            for b in ns:
                if a != b:
                    g[a].add(b)
        for a in ns:
            g[a].discard(v)
        del g[v]
        order.append(v)
    return width, order


def hist(sizes, label):
    import math
    print(f"\n  {label}: {len(sizes):,} blocks")
    bins = [(0, 500), (500, 1000), (1000, 2000), (2000, 3000),
            (3000, 4000), (4000, 5000), (5000, 10**9)]
    for lo, hi in bins:
        c = sum(1 for s in sizes if lo <= s < hi)
        if c:
            bar = '#' * max(1, int(40 * c / len(sizes)))
            tag = f"{lo:,}-{hi:,}" if hi < 10**9 else f">{lo:,}"
            print(f"    {tag:>14} : {c:6,} ({c/len(sizes):5.1%}) {bar}")
    ss = sorted(sizes)
    n = len(ss)
    print(f"    min={ss[0]:,}  p25={ss[n//4]:,}  median={ss[n//2]:,}  "
          f"p75={ss[3*n//4]:,}  max={ss[-1]:,}  total={sum(ss):,}")


def main():
    t0 = time.time()
    print("=" * 74)
    print("Q1 -- WHAT '80' WAS, AND WHAT ACTUALLY SETS SOLVE COST")
    print("=" * 74)
    print("""
  max |Q_ij| = 80  is a COEFFICIENT MAGNITUDE (7 bits of dynamic range).
  It says an annealer must distinguish 80 levels of coupler strength.
  It is NOT a number of coupled variables and implies NO 2^80 search.

  Exact solve cost is governed by the SEPARATOR of the block graph:
  dynamic programming over a tree decomposition costs about
        (#blocks) x 2^(bits in the largest separator).
  Below, separators are measured on the real block graph.
""")

    # ---- tiers ------------------------------------------------------------
    unk_b2, at_b2 = downstream(D2)
    unk_b1, at_b1 = downstream(DRIVERS)
    print(f"  tier [B''] : {len(at_b2)} atoms, {len(unk_b2)} unknown wires")
    print(f"  tier [B']  : {len(at_b1)} atoms, {len(unk_b1)} unknown wires")

    for label, atoms, unknown in [("[B''] branched", at_b2, unk_b2),
                                  ("[B'] constant-folded", at_b1, unk_b1)]:
        adj, holders = build_graph(atoms, unknown)
        deg = [len(v) for v in adj.values()]
        shared = sum(1 for x, s in holders.items() if len(s) > 1)
        twmd, _ = treewidth_ub(adj, 'min-degree')
        twmf, _ = treewidth_ub(adj, 'min-fill')
        tw = min(twmd, twmf)
        print(f"\n  --- {label} block graph ---")
        print(f"    blocks(nodes)      : {len(adj):,}")
        print(f"    shared wires       : {shared:,}  "
              f"(edges induced: {sum(deg)//2:,})")
        print(f"    block degree       : mean {sum(deg)/len(deg):.1f}, max {max(deg)}")
        print(f"    treewidth UB       : {tw} blocks "
              f"(min-degree {twmd}, min-fill {twmf})")
        # separator in BITS: a separator of t blocks carries the wires they share
        sep_bits = tw * W
        print(f"    separator in bits  : ~{tw} x {W} = ~{sep_bits:,} bits"
              f"  -> exact DP ~2^{sep_bits:,}")
        print(f"    => NOT 2^80; the honest exact-DP exponent here is "
              f"~2^{sep_bits:,}, i.e. infeasible.")

    # ---- Q2 size distribution ---------------------------------------------
    print("\n" + "=" * 74)
    print("Q2 -- SIZE DISTRIBUTION OF THE EMITTED QUBO BLOCKS")
    print("=" * 74)
    j = json.load(open(os.path.join(os.path.dirname(__file__), 'qubo_full.json')))
    print(f"  (from qubo_full.py census: q={j['q']}, w={j['w']})")
    # rebuild the packed size lists exactly as census does
    sys.path.insert(0, os.path.dirname(__file__))
    import qubo_full as QF
    mul_ref = QF.build_mul(W, j['q'])
    mul_size = mul_ref.nv - 3 * W
    lin_meas = {}
    for T in (2, 4, 8, 16):
        r = QF.build_lin(W, j['q'], [1] * T)
        lin_meas[T] = r.nv - T * W - W
    slope = (lin_meas[16] - lin_meas[2]) / 14.0
    base = lin_meas[2] - 2 * slope
    lin_cost = lambda T: int(base + slope * max(2, T))

    for label, atoms, unknown in [("[B''] branched", at_b2, unk_b2),
                                  ("[B'] constant-folded", at_b1, unk_b1)]:
        units = []
        for a in sorted(atoms):
            mm = ll = 0
            for mono in L.polys[a]:
                k = sum(1 for x in mono if x in unknown)
                if k >= 2: mm += 1
                elif k == 1: ll += 1
            units.extend([mul_size] * mm)
            units.append(lin_cost(ll))
        packed = QF.pack(units, 1000, 5000)
        hist(units, f"{label} -- RAW sub-blocks (before packing)")
        hist(packed, f"{label} -- PACKED QUBO blocks (target 1,000-5,000)")

    print(f"\n  {time.time()-t0:.1f}s")


if __name__ == '__main__':
    main()
