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


# ---------------------------------------------------- exact ceiling (replay)
# The encoding is a faithful sum-of-squares + Rosenberg-AND Hamiltonian, so for
# any (a,b,c) with c == a*b (mod p) there is exactly ONE zero-energy completion
# (every ancilla -- quotient, carries, adders, ANDs -- is a deterministic
# function computed by Q.witness).  Enumerating those completions over all valid
# triples therefore enumerates every ground state, in O(p * 2^s) replays instead
# of an exponential DFS.  (Cross-checked against verify.zero_states at p=13.)
def make(p, **kw):
    Q, A, B, C, _ = verify.make(p, **kw)
    return Q, A, B, C


def _triples(p, s, pin_a=None):
    """all ground-state triples: operands a,b are FREE s-bit words (range
    [0,2^s)), c is an s-bit word with c == a*b (mod p).  Each such triple has
    exactly one zero-energy completion (the encoding is a deterministic function
    of (a,b,c); verified by DFS at p=13), so iterating them enumerates every
    ground state.  pin_a fixes low bits of operand a: {bit:val}."""
    N = 1 << s
    for a in range(N):
        if pin_a and any(((a >> t) & 1) != v for t, v in pin_a.items()):
            continue
        for b in range(N):
            base = a * b % p
            c = base
            while c < N:
                yield a, b, c
                c += p


def _seen(Q, A, B, C, p, s, pin_a=None):
    seen0 = bytearray(Q.n)
    seen1 = bytearray(Q.n)
    cnt = 0
    for a, b, c in _triples(p, s, pin_a):
        inp = {}
        for t, v in enumerate(A.bits):
            inp[v] = (a >> t) & 1
        for t, v in enumerate(B.bits):
            inp[v] = (b >> t) & 1
        for t, v in enumerate(C.bits):
            inp[v] = (c >> t) & 1
        x, _ = Q.witness(inp, {'_a': a, '_b': b, '_c': c})
        for v in range(Q.n):
            if x[v]:
                seen1[v] = 1
            else:
                seen0[v] = 1
        cnt += 1
    return seen0, seen1, cnt


def ceiling(Q, A, B, C, p, s, pin_a=None):
    """variables constant over ALL ground states (optionally those with the
    given low bits of operand a pinned).  Returns (dict v->val, #states)."""
    seen0, seen1, cnt = _seen(Q, A, B, C, p, s, pin_a)
    fix = {v: (1 if seen1[v] else 0) for v in range(Q.n)
           if not (seen0[v] and seen1[v])}
    return fix, cnt


def verify_subset(fix, Q, A, B, C, p, s, pin_a=None, label=""):
    """gate: re-scan every ground state and assert each fixed var is constant."""
    for a, b, c in _triples(p, s, pin_a):
        inp = {}
        for t, v in enumerate(A.bits):
            inp[v] = (a >> t) & 1
        for t, v in enumerate(B.bits):
            inp[v] = (b >> t) & 1
        for t, v in enumerate(C.bits):
            inp[v] = (c >> t) & 1
        x, _ = Q.witness(inp, {'_a': a, '_b': b, '_c': c})
        for v, val in fix.items():
            assert x[v] == val, f"UNSOUND {label}: var {v} expected {val} got {x[v]} at a={a},b={b},c={c}"
    return True
