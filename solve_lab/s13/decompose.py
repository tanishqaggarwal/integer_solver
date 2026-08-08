#!/usr/bin/env python3
"""
Decompose the instance into minimally-coupled QUBO blocks of ~1,000-5,000
binary variables each.

Substrate: the RNS reformulation (S13_RNS_QUBO.md).  Work modulo a small prime
q with  w = ceil(log2 q)  bits per residue.  Each instance variable x_i is a
w-bit residue; each gate atom is a small penalty block over its variables'
residues plus multiply-auxiliaries.

Cut: the oriented circuit is a DAG (s9: 29,675 topo-ordered gates + checks).
Slicing the topological order into contiguous windows is the natural low-cut
partition for a DAG -- a variable produced in window k is consumed by gates that
come later in topo order, usually nearby, so few wires cross a window boundary.
Each window becomes one QUBO block; variables consumed from other windows are the
INTERFACE (the coupling), held fixed during a block solve (block coordinate
descent / ADMM over the blocks).

Concrete QUBO cost model (per block b), stated so "1,000-5,000 variables" is
well-defined:
    cost(b) = w * |owned vars in b|                      # residue registers
            + sum over gate atoms g in b of AUX(deg g)   # multiply auxiliaries
    AUX(1) = 0                      # linear atom: no product
    AUX(2) = w*w + 3*w              # one w x w schoolbook multiply + mod-q reduce
    AUX(4) = 3*(w*w + 3*w)          # degree-4 atom ~ 3 multiplies
Interface variables are inputs, not free QUBO variables, so they are NOT counted
in cost(b); they are the coupling we minimize and report.

Usage:  python3 decompose.py [w] [target_lo] [target_hi]
        (defaults w=8, target 1000..5000 binary vars/block)
Writes:  s13/blocks_w{w}.json  (block manifest) and prints coupling metrics.
"""
import os, sys, json, time
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 's9', 'eff'))
import lib as L

W        = int(sys.argv[1]) if len(sys.argv) > 1 else 8
TGT_LO   = int(sys.argv[2]) if len(sys.argv) > 2 else 1000
TGT_HI   = int(sys.argv[3]) if len(sys.argv) > 3 else 5000
OUT      = os.path.join(os.path.dirname(__file__), f'blocks_w{W}.json')

NA  = len(L.avars)                       # atoms
topo = list(L.topo)                      # gate atoms in dependency order
topo_pos = {a: i for i, a in enumerate(topo)}


def degree(atom):
    d = 0
    for mono in L.polys[atom]:
        d = max(d, len(mono))
    return d

AUX = {0: 0, 1: 0, 2: W * W + 3 * W, 4: 3 * (W * W + 3 * W)}
def aux(atom):
    d = degree(atom)
    return AUX.get(d, 3 * (W * W + 3 * W))


def build_blocks(nblocks):
    """Assign every atom and every variable to one of nblocks blocks."""
    T = len(topo)
    atom_block = {}
    # 1) gate atoms in topo order -> contiguous windows
    for i, a in enumerate(topo):
        atom_block[a] = min(nblocks - 1, i * nblocks // T)
    # 2) owner block of each variable produced by a topo gate
    owner = {}
    for a in topo:
        oc = L.atom_out.get(a)
        if oc is not None:
            owner[oc[1]] = atom_block[a]
    # 3) remaining atoms (checks + cyclic gates): block = max owner over their
    #    already-placed input vars (so a check sits with its latest input);
    #    variables with no owner yet are resolved after.
    unplaced = [a for a in range(NA) if a not in atom_block]
    # first, free inputs / cyclic vars: owner = min block of any topo gate using them
    var_use_block = {}
    for a in topo:
        b = atom_block[a]
        for v in L.avars[a]:
            if v not in owner:
                if v not in var_use_block or b < var_use_block[v]:
                    var_use_block[v] = b
    for v, b in var_use_block.items():
        owner.setdefault(v, b)
    # place remaining atoms
    for a in unplaced:
        bs = [owner[v] for v in L.avars[a] if v in owner]
        atom_block[a] = max(bs) if bs else 0
    # any var still without owner -> block of an atom that references it
    for a in range(NA):
        b = atom_block[a]
        for v in L.avars[a]:
            owner.setdefault(v, b)
    return atom_block, owner


def measure(atom_block, owner, nblocks):
    owned = [0] * nblocks
    cost = [0] * nblocks
    gates = [0] * nblocks
    checks = [0] * nblocks
    for v, b in owner.items():
        owned[b] += 1
        cost[b] += W
    for a in range(NA):
        b = atom_block[a]
        if a in topo_pos or L.atom_out.get(a) is not None:
            gates[b] += 1
            cost[b] += aux(a)
        else:
            checks[b] += 1
    # coupling: for each atom, count input vars owned elsewhere
    interface_incidences = 0
    boundary_vars = set()
    fanout = {}                       # var -> set of blocks that consume it
    for a in range(NA):
        b = atom_block[a]
        for v in L.avars[a]:
            fanout.setdefault(v, set()).add(b)
            if owner[v] != b:
                interface_incidences += 1
                boundary_vars.add(v)
    spanning = sum(1 for v, bs in fanout.items() if len(bs) > 1)
    return dict(owned=owned, cost=cost, gates=gates, checks=checks,
                interface_incidences=interface_incidences,
                boundary_vars=len(boundary_vars), spanning_vars=spanning,
                fanout=fanout)


def main():
    t0 = time.time()
    print(f"[decompose] w={W} bits/residue; target {TGT_LO}-{TGT_HI} binary vars/block")
    # pick nblocks so the MEAN block cost lands mid-target; then verify spread
    total_cost = sum(W for _ in owner_all()) + sum(aux(a) for a in range(NA)
                                                   if L.atom_out.get(a) is not None)
    mid = (TGT_LO + TGT_HI) // 2
    nblocks = max(1, round(total_cost / mid))
    print(f"[decompose] total QUBO cost ~= {total_cost:,} binary vars "
          f"-> {nblocks} blocks of ~{total_cost//nblocks} each")

    ab, owner = build_blocks(nblocks)
    m = measure(ab, owner, nblocks)
    costs = m['cost']
    inrange = sum(1 for c in costs if TGT_LO <= c <= TGT_HI)
    print(f"\n[decompose] {nblocks} blocks built in {time.time()-t0:.1f}s")
    print(f"  block QUBO size  : min={min(costs)}  mean={sum(costs)//len(costs)}  max={max(costs)}")
    print(f"  in target range  : {inrange}/{nblocks} blocks within [{TGT_LO},{TGT_HI}]")
    print(f"  gates/block      : mean={sum(m['gates'])//nblocks}")
    print(f"  COUPLING:")
    print(f"    interface incidences (input wires crossing a block edge): {m['interface_incidences']:,}")
    print(f"    unique boundary variables                              : {m['boundary_vars']:,}")
    print(f"    variables spanning >1 block                            : {m['spanning_vars']:,}")
    tot_var_incidence = sum(len(L.avars[a]) for a in range(NA))
    print(f"    interface fraction (incidences / all var-incidences)   : "
          f"{m['interface_incidences']/tot_var_incidence:.3%}")
    # hubs: variables consumed by the most blocks -- the irreducible coupling
    fo = m['fanout']
    hubs = sorted(fo.items(), key=lambda kv: -len(kv[1]))[:12]
    print(f"  worst coupling hubs (var: #blocks consuming it):")
    for v, bs in hubs:
        print(f"    x_{v}: {len(bs)} blocks")

    manifest = {
        'w_bits': W, 'nblocks': nblocks,
        'target': [TGT_LO, TGT_HI],
        'block_cost': costs, 'block_owned': m['owned'],
        'block_gates': m['gates'], 'block_checks': m['checks'],
        'interface_incidences': m['interface_incidences'],
        'boundary_vars': m['boundary_vars'],
        'spanning_vars': m['spanning_vars'],
        'hubs': [[v, len(bs)] for v, bs in hubs],
    }
    json.dump(manifest, open(OUT, 'w'))
    print(f"\n[decompose] manifest -> {OUT}  ({time.time()-t0:.1f}s)")


def owner_all():
    # helper: count of variables that appear anywhere (for cost pre-estimate)
    seen = set()
    for a in range(NA):
        seen |= L.avars[a]
    return seen


if __name__ == '__main__':
    main()
