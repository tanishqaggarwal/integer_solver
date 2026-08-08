#!/usr/bin/env python3
"""persist.py -- persistency: what a presolver can PROVABLY fix.

Soundness contract (the gate): a variable may be reported fixed to v only if v is
its value in EVERY ground state.  Every method here is checked against the
exhaustively enumerated ground-state set (verify.zero_states) on small instances.

Because every penalty is a perfect square / Rosenberg AND penalty, the ground
states are exactly E=0, i.e. the SOLUTIONS of the constraint system.  This is a
pure satisfaction problem, so:

  * WEAK persistency (a value that occurs in SOME global minimum) is vacuous --
    for any extendable variable both values occur in some solution, so it fixes
    nothing.  Classic QPBO 'autarky' yields only weak persistency; here it is
    empty.  We verify this.

  * STRONG persistency (value constant over ALL global minima) is the only thing
    that reduces the problem.  It is EXACTLY {v : v constant over all ground
    states} = the enumerated ceiling.  Roof duality / QPBO strong persistency is
    a SUBSET of this ceiling (a theorem), so the ceiling upper-bounds every
    presolver.  Probing+propagation is our sound constructive lower bound.
"""
import os, sys
from collections import defaultdict, deque

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import common  # noqa: E402
sys.path.insert(0, common.SQUEEZE)
from presolve import propagate, gates_of  # noqa: E402
import verify  # noqa: E402


# ---------------------------------------------------- constructive presolvers
def prop_fix(Q, fixed=None):
    """unit/interval propagation to a fixpoint. Sound. Returns dict v->val."""
    return propagate(Q.squares, gates_of(Q), Q.n, fixed=dict(fixed or {}))


def probe_fix(Q, fixed=None, candidates=None, rounds=3):
    """probing: fix each candidate var to 0 and to 1, propagate; if one branch is
    infeasible fix the other, if BOTH feasible intersect the two derived domains
    and keep every variable forced to the same value in both.  Sound (both are
    necessary-consequence closures). Iterated to a fixpoint.  Returns dict."""
    gates = gates_of(Q)
    dom = dict(prop_fix(Q, fixed))
    if candidates is None:
        candidates = list(range(Q.n))
    for _ in range(rounds):
        moved = False
        for v in candidates:
            if v in dom:
                continue
            try:
                d0 = propagate(Q.squares, gates, Q.n, fixed={**dom, v: 0})
            except ValueError:
                d0 = None
            try:
                d1 = propagate(Q.squares, gates, Q.n, fixed={**dom, v: 1})
            except ValueError:
                d1 = None
            if d0 is None and d1 is None:
                raise ValueError("infeasible")
            if d0 is None:
                dom = d1
                moved = True
                continue
            if d1 is None:
                dom = d0
                moved = True
                continue
            # both feasible: keep agreements
            for u, val in d0.items():
                if u not in dom and d1.get(u) == val:
                    dom[u] = val
                    moved = True
        if not moved:
            break
    return dom


# ---------------------------------------------------- exact ceiling (enumerate)
def ground_states(p, **kw):
    """every (full variable assignment) ground state, plus Q,A,B,C."""
    Q, A, B, C, _ = verify.make(p, **kw)
    xs = verify.zero_states(Q)
    return Q, A, B, C, xs


def ceiling(Q, xs, restrict=None):
    """variables constant over the ground-state list xs (optionally the subset
    consistent with the partial assignment `restrict`). Returns dict v->val."""
    if restrict:
        xs = [x for x in xs if all(x[v] == val for v, val in restrict.items())]
    if not xs:
        return {}, 0
    n = Q.n
    seen0 = [False] * n
    seen1 = [False] * n
    for x in xs:
        for v in range(n):
            if x[v]:
                seen1[v] = True
            else:
                seen0[v] = True
    fix = {v: (1 if seen1[v] else 0) for v in range(n)
           if not (seen0[v] and seen1[v])}
    return fix, len(xs)


def verify_subset(fix, Q, xs, label=""):
    """gate: every fixed var must hold that value in EVERY ground state."""
    bad = []
    for v, val in fix.items():
        for x in xs:
            if x[v] != val:
                bad.append((v, val))
                break
    assert not bad, f"UNSOUND {label}: {bad[:10]}"
    return True
