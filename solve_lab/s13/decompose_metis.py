#!/usr/bin/env python3
"""
Minimally-coupled QUBO decomposition via METIS multilevel k-way partitioning.

Model (RNS substrate, S13_RNS_QUBO.md): work mod a small prime q, w=ceil(log2 q)
bits per residue.  Build a graph whose NODES are the instance variables and whose
EDGES join variables that co-occur in an atom (the gate/check that couples them).

  * Node weight = w * (residue bits) + the multiply-auxiliary bits of the atom
    that DEFINES the variable, so a balanced partition balances QUBO binary-var
    count, not raw variable count.
  * Edge weight = number of atoms coupling the two variables.

METIS minimizes total cut edge weight for a balanced k-way partition -> each part
is a QUBO block; a cut edge is an inter-block coupling (an interface variable that
one block produces and another consumes, held fixed during a block solve / ADMM).

We choose k so the mean block lands mid-target and report the achieved size spread
and coupling.  Compared with naive topological slicing (decompose.py) the interface
fraction drops from ~66% to the value printed here.

Usage:  python3 decompose_metis.py [w] [target_lo] [target_hi]
Writes: s13/blocks_metis_w{w}.json
"""
import os, sys, json, time
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 's9', 'eff'))
import lib as L
import pymetis

W      = int(sys.argv[1]) if len(sys.argv) > 1 else 8
TGT_LO = int(sys.argv[2]) if len(sys.argv) > 2 else 1000
TGT_HI = int(sys.argv[3]) if len(sys.argv) > 3 else 5000
OUT    = os.path.join(os.path.dirname(__file__), f'blocks_metis_w{W}.json')
NA = len(L.avars)

AUX = {0: 0, 1: 0, 2: W * W + 3 * W, 4: 3 * (W * W + 3 * W)}
def degree(atom):
    return max((len(m) for m in L.polys[atom]), default=0)
def aux(atom):
    return AUX.get(degree(atom), 3 * (W * W + 3 * W))


def main():
    t0 = time.time()
    print(f"[metis] w={W} bits/residue; target {TGT_LO}-{TGT_HI} binary vars/block")

    # ---- node weights: residue bits + aux of the defining atom ----------------
    allvars = sorted(set().union(*L.avars))
    idx = {v: i for i, v in enumerate(allvars)}
    n = len(allvars)
    nodew = [W] * n                       # every residue costs w bits
    for a in range(NA):
        oc = L.atom_out.get(a)
        if oc is not None and oc[1] in idx:
            nodew[idx[oc[1]]] += aux(a)   # its producing gate's multiply aux

    # ---- adjacency: variables co-occurring in an atom (weighted) --------------
    from collections import defaultdict
    ew = defaultdict(int)
    for a in range(NA):
        vs = [idx[v] for v in L.avars[a]]
        for i in range(len(vs)):
            for j in range(i + 1, len(vs)):
                x, y = vs[i], vs[j]
                if x > y: x, y = y, x
                ew[(x, y)] += 1
    adj = [[] for _ in range(n)]
    adjw = [[] for _ in range(n)]
    for (x, y), c in ew.items():
        adj[x].append(y); adjw[x].append(c)
        adj[y].append(x); adjw[y].append(c)
    print(f"[metis] graph: {n:,} nodes  {len(ew):,} edges  (built {time.time()-t0:.1f}s)")

    total_cost = sum(nodew)
    k = max(2, round(total_cost / ((TGT_LO + TGT_HI) // 2)))
    print(f"[metis] total QUBO cost ~= {total_cost:,} bits -> k={k} blocks "
          f"of ~{total_cost//k} each")

    # ---- METIS k-way, balancing node weight (QUBO cost) -----------------------
    xadj = [0]; adjncy = []; eweights = []
    for i in range(n):
        adjncy.extend(adj[i]); eweights.extend(adjw[i]); xadj.append(len(adjncy))
    t1 = time.time()
    cutcount, membership = pymetis.part_graph(
        k, xadj=xadj, adjncy=adjncy, vweights=nodew, eweights=eweights)
    membership = list(membership)
    print(f"[metis] partitioned in {time.time()-t1:.1f}s; METIS edge-cut={cutcount:,}")

    # ---- measure --------------------------------------------------------------
    cost = [0] * k; owned = [0] * k
    for i in range(n):
        b = membership[i]; owned[b] += 1; cost[b] += nodew[i]
    # coupling: variable spans >1 block if some atom puts it with a different owner
    fanout = defaultdict(set)
    for a in range(NA):
        for v in L.avars[a]:
            fanout[idx[v]].add(membership[idx[v]])   # own block always
    # a var is an interface var iff it is CONSUMED in a block != its own
    interface_inc = 0; boundary = set()
    for a in range(NA):
        bs = set(membership[idx[v]] for v in L.avars[a])
        if len(bs) > 1:
            # this atom is cut; every var whose owner != this atom's home couples
            home = membership[idx[L.atom_out[a][1]]] if L.atom_out.get(a) and L.atom_out[a][1] in idx else min(bs)
            for v in L.avars[a]:
                if membership[idx[v]] != home:
                    interface_inc += 1; boundary.add(v)
    costs_sorted = sorted(cost)
    inrange = sum(1 for c in cost if TGT_LO <= c <= TGT_HI)
    tot_inc = sum(len(L.avars[a]) for a in range(NA))
    span = sum(1 for v, bs in fanout.items() if len(bs) > 1)  # placeholder (own only)
    # true spanning: count vars consumed across a cut
    consumed_blocks = defaultdict(set)
    for a in range(NA):
        b_atom = membership[idx[L.atom_out[a][1]]] if L.atom_out.get(a) and L.atom_out[a][1] in idx else None
        for v in L.avars[a]:
            consumed_blocks[v].add(membership[idx[v]])
            if b_atom is not None:
                consumed_blocks[v].add(b_atom)
    spanning = sum(1 for v, bs in consumed_blocks.items() if len(bs) > 1)

    print(f"\n[metis] {k} blocks")
    print(f"  block QUBO size : min={min(cost)}  p10={costs_sorted[k//10]}  "
          f"mean={sum(cost)//k}  p90={costs_sorted[9*k//10]}  max={max(cost)}")
    print(f"  in target range : {inrange}/{k} blocks within [{TGT_LO},{TGT_HI}]")
    print(f"  COUPLING:")
    print(f"    METIS edge-cut (weighted)        : {cutcount:,}")
    print(f"    interface incidences             : {interface_inc:,}")
    print(f"    unique boundary variables        : {len(boundary):,} of {n:,} "
          f"({len(boundary)/n:.1%})")
    print(f"    variables spanning >1 block      : {spanning:,} ({spanning/n:.1%})")
    print(f"    interface fraction               : {interface_inc/tot_inc:.2%}")

    hubs = sorted(consumed_blocks.items(), key=lambda kv: -len(kv[1]))[:12]
    print(f"  worst coupling hubs (var: #blocks touching it):")
    for v, bs in hubs:
        print(f"    x_{v}: {len(bs)}")

    json.dump({
        'w_bits': W, 'k': k, 'target': [TGT_LO, TGT_HI],
        'metis_cut': cutcount, 'block_cost': cost, 'block_owned': owned,
        'interface_incidences': interface_inc,
        'boundary_vars': len(boundary), 'spanning_vars': spanning,
        'membership': membership, 'allvars': allvars,
        'hubs': [[int(v), len(bs)] for v, bs in hubs],
    }, open(OUT, 'w'))
    print(f"\n[metis] manifest -> {OUT}  ({time.time()-t0:.1f}s)")


if __name__ == '__main__':
    main()
