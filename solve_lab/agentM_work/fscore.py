"""Fast exact equation scoring.

E2.eqfails rescans all 39,033 equations for every candidate.  Only equations that
touch a currently-nonzero atom can differ from the constant-only baseline, so score
in O(#touched equations) instead.  Exact, not an approximation.
"""
import sys, os, collections
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentE_work')
sys.set_int_max_str_digits(20_000_000)
import harness as H

NEQ = len(H.eqt)
_const = [0] * NEQ
_cmap = [None] * NEQ
atom2eq = collections.defaultdict(list)
for e, (issq, outer, terms) in enumerate(H.eqt):
    c0 = 0
    cm = {}
    for c, a in terms:
        if a < 0:
            c0 += c
        else:
            cm[a] = cm.get(a, 0) + c
    _const[e] = c0
    _cmap[e] = cm
    for a in cm:
        atom2eq[a].append(e)

# equations that fail when NO atom is nonzero
_BASEFAIL = set(e for e in range(NEQ) if _const[e] != 0)


def fails(bad):
    """Exact set of failing equation indices for a bad-atom dict."""
    cand = set()
    for a in bad:
        cand.update(atom2eq[a])
    out = set(e for e in _BASEFAIL if e not in cand)
    for e in cand:
        s = _const[e]
        cm = _cmap[e]
        for a, v in bad.items():
            c = cm.get(a)
            if c:
                s += c * v
        if s:
            out.add(e)
    return out


def nfail(bad):
    cand = set()
    for a in bad:
        cand.update(atom2eq[a])
    n = len(_BASEFAIL - cand)
    for e in cand:
        s = _const[e]
        cm = _cmap[e]
        for a, v in bad.items():
            c = cm.get(a)
            if c:
                s += c * v
        if s:
            n += 1
    return n


def score(bad):
    return NEQ - nfail(bad)
