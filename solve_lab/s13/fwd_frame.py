#!/usr/bin/env python3
"""
Block-preserving forward evaluator (the 'frame' the lab's RESUME warns about).

The 39,026 witness is OFF-MANIFOLD: two gate atoms are deliberately nonzero
(a35761 -> x31864 and a22229 -> x7068).  A naive forward pass "repairs" them and
silently changes 109 downstream wires, which is why the first linearity fit was
computed on garbage.

Correct rule: a wire is recomputed ONLY if its defining gate is satisfied at the
base state.  Outputs of broken gates are FROZEN at their stored values and act
as extra inputs.  Everything else propagates in topological order.

evaluate(v, overrides, unknown, touched) -> dict wire -> value
    v         : the base assignment (list)
    overrides : {wire: new value} for the free inputs being moved
"""
import os, sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 's9', 'eff'))
import lib as L


def broken_gates(v, touched=None):
    """Gate atoms that are NONZERO at v -- their outputs must not be recomputed."""
    av = L.all_atom_values(v)
    out = {}
    src = touched if touched is not None else range(len(L.avars))
    for a in src:
        oc = L.atom_out.get(a)
        if oc is not None and av[a] != 0:
            out[oc[1]] = a
    return out


def eval_atom_without(a, y, val, v):
    """Evaluate atom a with the monomials containing y dropped."""
    tot = 0
    for mono, c in L.polys[a].items():
        if y in mono:
            continue
        t = c
        for x in mono:
            t *= val[x] if x in val else v[x]
        tot += t
    return tot


def evaluate(v, overrides, unknown, touched, frozen=None):
    """
    Propagate `overrides` through the cone, preserving broken gates.
    Returns {wire: value} for every wire in `unknown`.
    """
    if frozen is None:
        frozen = broken_gates(v, touched)

    val = {}
    for x in unknown:
        val[x] = v[x]                    # start from the state
    val.update(overrides)

    # which wires are genuinely recomputable
    recompute = {}
    for a in touched:
        oc = L.atom_out.get(a)
        if oc is None:
            continue
        y = oc[1]
        if y in unknown and y not in overrides and y not in frozen:
            recompute[y] = a

    # topological order restricted to the cone
    order = [y for y in L.topo if y in recompute] if hasattr(L, 'topo') else []
    seen = set(order)
    order += [y for y in recompute if y not in seen]

    # L.topo lists ATOMS in some builds and VARS in others; handle both by
    # falling back to iterate-to-fixpoint over the recompute set.
    for _ in range(len(recompute) + 3):
        changed = False
        for y, a in recompute.items():
            coeff = L.atom_out[a][0]
            tot = eval_atom_without(a, y, val, v)
            newv = -tot // coeff if coeff in (1, -1) else None
            if newv is None:
                continue
            if val.get(y) != newv:
                val[y] = newv; changed = True
        if not changed:
            break
    return val


def check_value(a, v, val):
    tot = 0
    for mono, c in L.polys[a].items():
        t = c
        for x in mono:
            t *= val[x] if x in val else v[x]
        tot += t
    return tot


def self_test(v, unknown, touched):
    """Identity check: no overrides must reproduce the state exactly."""
    frozen = broken_gates(v, touched)
    val = evaluate(v, {}, unknown, touched, frozen)
    bad = [x for x in unknown if val[x] != v[x]]
    return frozen, bad
