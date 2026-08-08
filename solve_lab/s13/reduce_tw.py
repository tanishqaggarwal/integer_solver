#!/usr/bin/env python3
"""
ATTACK THE COUPLING: can the treewidth between blocks actually be reduced?

The earlier number (blockgraph.py) was measured on the ATOM graph and then
multiplied by 16 bits/wire.  That is the wrong object twice over:

  * the constraint graph for a tree decomposition is over WIRES, not atoms;
  * more importantly, a wire that is DETERMINED by a gate is not a degree of
    freedom.  It can be SUBSTITUTED AWAY.  Carrying it in a DP state is pure
    waste -- the correct state is only the FREE inputs that are still live.

This script does the reduction properly, in four stages, measuring each:

  S1  wire-level constraint graph + treewidth (the honest baseline)
  S2  eliminate determined wires by substitution -> effective graph over FREE
      inputs only.  This is the reduction the earlier analysis missed.
  S3  SPLIT test: does each constraint depend on x9118, on x8731, or on BOTH?
      If few constraints touch both, the problem separates into two independent
      256-bit problems and the coupling collapses.
  S4  LINEARITY test: is each constraint linear in the unknowns mod p?  A linear
      system is solved by elimination -- treewidth becomes irrelevant.

Usage: python3 reduce_tw.py
"""
import os, sys, time, json
from collections import defaultdict, deque

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L

LAB = os.path.join(HERE, '..')
P = 2**256 - 2**32 - 977
NA = len(L.avars)
D2 = [8731, 9118]


def treewidth_ub(adj, method='min-fill'):
    g = {v: set(ns) for v, ns in adj.items()}
    width = 0
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
    return width


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


# ---------------------------------------------------------------- evaluator --
def local_forward(v, free_vals, unknown, touched):
    """
    Recompute every DETERMINED unknown wire from the free inputs, in dependency
    order, using its own defining gate (output coefficient is +-1).
    Returns a dict wire -> value (mod-free, exact integers).
    """
    val = dict(free_vals)
    # topological order restricted to the cone
    defined = {}
    for a in touched:
        oc = L.atom_out.get(a)
        if oc is not None and oc[1] in unknown and oc[1] not in free_vals:
            defined[oc[1]] = a
    # iterate to a fixpoint (cone is small and acyclic in practice)
    pending = set(defined)
    for _ in range(len(defined) + 5):
        progressed = False
        for y in list(pending):
            a = defined[y]
            coeff = L.atom_out[a][0]
            # evaluate the atom with y omitted
            tot = 0
            ready = True
            for mono, c in L.polys[a].items():
                if y in mono:
                    continue
                t = c
                for x in mono:
                    xv = val.get(x, v[x] if x not in unknown else None)
                    if xv is None:
                        ready = False; break
                    t *= xv
                if not ready:
                    break
                tot += t
            if not ready:
                continue
            # atom == 0  =>  coeff*y + tot == 0  =>  y = -tot/coeff
            if tot % coeff != 0:
                val[y] = -tot // coeff
            else:
                val[y] = -tot // coeff
            pending.discard(y); progressed = True
        if not pending or not progressed:
            break
    return val


def check_value(a, v, val, unknown):
    tot = 0
    for mono, c in L.polys[a].items():
        t = c
        for x in mono:
            t *= val.get(x, v[x])
        tot += t
    return tot


def main():
    t0 = time.time()
    v = L.load(os.path.join(LAB, 'best', 'new_instance_partial_39026.json'))
    unknown, touched = downstream(D2)
    checks = sorted(a for a in touched if L.atom_out.get(a) is None)
    print(f"cone: {len(unknown)} unknown wires, {len(touched)} atoms, "
          f"{len(checks)} checks\n")

    # ---------------- S1 : wire-level constraint graph -----------------------
    print("=" * 74)
    print("S1 -- the HONEST baseline: treewidth of the WIRE constraint graph")
    print("=" * 74)
    adj = defaultdict(set)
    for a in touched:
        vs = [x for x in L.avars[a] if x in unknown]
        for i in range(len(vs)):
            for j in range(i + 1, len(vs)):
                adj[vs[i]].add(vs[j]); adj[vs[j]].add(vs[i])
    for x in unknown:
        adj.setdefault(x, set())
    tw1 = min(treewidth_ub(adj, 'min-fill'), treewidth_ub(adj, 'min-degree'))
    print(f"  nodes (wires)      : {len(adj)}")
    print(f"  treewidth UB       : {tw1} wires")
    print(f"  naive DP state     : {tw1} x 256 bits = {tw1*256:,} bits")
    print(f"  (the earlier atom-graph number was the wrong object)")

    # ---------------- S2 : eliminate determined wires ------------------------
    print("\n" + "=" * 74)
    print("S2 -- REDUCTION: substitute away every DETERMINED wire")
    print("=" * 74)
    free_all = {t for t in range(L.NVARS) if t not in L.definer}
    free_in_cone = sorted(x for x in unknown if x in free_all)
    derived = [x for x in unknown if x not in free_all]
    print(f"  wires in cone        : {len(unknown)}")
    print(f"  DETERMINED by a gate : {len(derived)}  -> eliminated, cost 0")
    print(f"  genuinely FREE       : {len(free_in_cone)} "
          f"{['x%d' % x for x in free_in_cone]}")
    print(f"  => after substitution every constraint is a function of only")
    print(f"     these {len(free_in_cone)} unknowns.  Effective treewidth <= "
          f"{len(free_in_cone)}, i.e. <= {len(free_in_cone)*256} bits,")
    print(f"     NOT {tw1*256:,}.  The determined wires were never state.")

    # ---------------- S3 : does the problem SPLIT? ---------------------------
    print("\n" + "=" * 74)
    print("S3 -- SPLIT TEST: does each constraint touch one unknown or both?")
    print("=" * 74)
    # structural support: which free inputs each wire depends on
    supp = {x: ({x} if x in free_all else set()) for x in unknown}
    for _ in range(len(unknown) + 2):
        changed = False
        for a in touched:
            oc = L.atom_out.get(a)
            if oc is None or oc[1] not in unknown:
                continue
            y = oc[1]
            if y in free_all:
                continue
            s = set()
            for x in L.avars[a]:
                if x == y:
                    continue
                s |= supp.get(x, set())
            if not s <= supp[y]:
                supp[y] |= s; changed = True
        if not changed:
            break
    only_a = only_b = both = neither = 0
    both_list = []
    for a in checks:
        s = set()
        for x in L.avars[a]:
            s |= supp.get(x, set())
        ha, hb = (8731 in s), (9118 in s)
        if ha and hb:
            both += 1; both_list.append(a)
        elif ha:
            only_a += 1
        elif hb:
            only_b += 1
        else:
            neither += 1
    print(f"  checks depending on x8731 ONLY : {only_a}")
    print(f"  checks depending on x9118 ONLY : {only_b}")
    print(f"  checks depending on BOTH       : {both}  {['a%d'%a for a in both_list[:12]]}")
    print(f"  checks depending on NEITHER    : {neither}")
    if both == 0:
        print(f"  ==> THE PROBLEM SPLITS COMPLETELY into two independent")
        print(f"      256-bit subproblems.  Coupling between them is ZERO.")
    else:
        print(f"  ==> {both} constraint(s) couple the two unknowns; everything")
        print(f"      else is separable.  The coupling is exactly those {both}.")

    # ---------------- S4 : linearity in the unknowns -------------------------
    print("\n" + "=" * 74)
    print("S4 -- LINEARITY: is each constraint linear in the unknowns (mod p)?")
    print("=" * 74)
    print("  (a linear system is solved by elimination; treewidth is then moot)")
    base = {x: v[x] for x in D2}
    results = defaultdict(int)
    lin_checks = []
    for a in checks[:60]:
        vals = []
        okrow = True
        for t in (0, 1, 2):
            fv = dict(base)
            fv[9118] = base[9118] + t
            try:
                val = local_forward(v, fv, unknown, touched)
                vals.append(check_value(a, v, val, unknown) % P)
            except Exception:
                okrow = False; break
        if not okrow or len(vals) < 3:
            results['error'] += 1; continue
        second = (vals[2] - 2 * vals[1] + vals[0]) % P
        if second == 0:
            results['linear in x9118'] += 1
            lin_checks.append(a)
        else:
            results['nonlinear in x9118'] += 1
    for k, n in sorted(results.items()):
        print(f"    {k:24s}: {n}")
    print(f"  => if the constraints are linear mod p, the residual is a LINEAR")
    print(f"     SYSTEM over GF(p) in {len(free_in_cone)} unknowns -- Gaussian")
    print(f"     elimination, no search, no annealer, treewidth irrelevant.")

    print(f"\n  {time.time()-t0:.1f}s")


if __name__ == '__main__':
    main()
