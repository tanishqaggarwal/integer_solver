"""Corrected engine: demote the atoms that E's orientation wrongly uses as definers.

E's _bootstrap picks, for each derived var u, an atom i that it then SOLVES TO ZERO.
Any atom used as a definer is therefore identically satisfied in every state E's
forward can produce.  The 39,026 deliverable deliberately leaves 8 atoms nonzero,
5 of which are definers in E's orientation -- so E's reachable space excludes the
deliverable entirely.  Fix: demote those atoms, promoting their vars to free inputs.
"""
import sys, os, math, json, pickle, collections
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentE_work')
sys.set_int_max_str_digits(20_000_000)
import harness as H

NV = H.NV
atoms = H.atoms; acodes = H.acodes; avars = H.avars; occ = H.occ

# ---- the demotion set -------------------------------------------------------
# atoms E uses as definers but which the best-known point needs nonzero
DEMOTE_ATOMS = [23616, 23617, 36659, 36663, 36664]

_atom2var = {}
for _u in H.SEQ:
    _i, _k = H.definer[_u]
    _atom2var[_i] = _u

PIN = sorted(_atom2var[a] for a in DEMOTE_ATOMS if a in _atom2var)  # vars promoted to free
PINSET = set(PIN)

definer = list(H.definer)
for u in PIN:
    definer[u] = None
SEQ = [u for u in H.SEQ if u not in PINSET]
FREE = sorted(set(H.FREE) | PINSET | {u for u in range(NV) if definer[u] is None})

SOLVE = [(u, definer[u][0], definer[u][1][0]) for u in SEQ]


def _solvevar(v, ns, u, i, kd):
    c = acodes[i]
    v[u] = 0; c0 = eval(c, ns)
    v[u] = 1; c1 = eval(c, ns)
    if kd == 'lin':
        sl = c1 - c0
        v[u] = -c0 // sl if sl and c0 % sl == 0 else 0
    else:
        v[u] = 2; c2 = eval(c, ns)
        A2 = c2 - 2 * c1 + c0; A = A2 // 2; B = c1 - c0 - A; C = c0
        disc = B * B - 4 * A * C
        if disc < 0 or A == 0:
            v[u] = 0; return
        r = math.isqrt(disc)
        if r * r != disc:
            v[u] = 0; return
        rts = {(-B + s) // (2 * A) for s in (r, -r) if (-B + s) % (2 * A) == 0}
        v[u] = rts.pop() if len(rts) == 1 else 0


def forward(seed):
    v = [0] * NV
    for k, val in seed.items():
        v[k] = val
    ns = {'v': v, '__builtins__': {}}
    for u, i, kd in SOLVE:
        _solvevar(v, ns, u, i, kd)
    return v


def badatoms(v):
    ns = {'v': v, '__builtins__': {}}
    out = {}
    for i in range(len(atoms)):
        r = eval(acodes[i], ns)
        if r:
            out[i] = r
    return out


def eqfails(av):
    fails = []
    for e, (issq, outer, terms) in enumerate(H.eqt):
        s = 0
        for c, a in terms:
            if a < 0:
                s += c
            elif a in av:
                s += c * av[a]
        if s:
            fails.append(e)
    return fails


def score(seed):
    v = forward(seed)
    av = badatoms(v)
    return len(eqfails(av)), av, v


def seed_of(v):
    """Extract the free-input seed (incl. pinned vars) from a full vector."""
    return {f: v[f] for f in FREE if v[f] != 0}
