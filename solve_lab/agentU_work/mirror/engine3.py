"""Engine with a CONFIGURABLE demotion set, for pricing alternative defect placements.

engine2 demotes a fixed 5 atoms (the ones the 39,026 point needs nonzero).  To search
*placements* rather than repairs we need to demote arbitrary atom sets: demoting atom a
frees the variable it defines, which turns a into a quantity we can move.

Key property that makes the search safe: demoting an atom and seeding its variable with
its CURRENT value leaves the whole state bit-identical, so demotion is score-neutral and
purely adds a degree of freedom.  Verified by `validate()`.
"""
import sys, os, math, json, collections
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentU_work/mirror')
sys.set_int_max_str_digits(20_000_000)
import harness as H

NV = H.NV
atoms = H.atoms; acodes = H.acodes; avars = H.avars; occ = H.occ

ATOM2VAR = {}
for _u in H.SEQ:
    _i, _k = H.definer[_u]
    ATOM2VAR[_i] = _u

BASE_DEMOTE = [23616, 23617, 36659, 36663, 36664]      # engine2's set


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


class Eng:
    def __init__(self, demote):
        self.demote = sorted(set(demote))
        self.pin = sorted(ATOM2VAR[a] for a in self.demote if a in ATOM2VAR)
        pinset = set(self.pin)
        self.definer = list(H.definer)
        for u in self.pin:
            self.definer[u] = None
        self.SEQ = [u for u in H.SEQ if u not in pinset]
        self.FREE = sorted(set(H.FREE) | pinset |
                           {u for u in range(NV) if self.definer[u] is None})
        self.SOLVE = [(u, self.definer[u][0], self.definer[u][1][0]) for u in self.SEQ]
        self._pos = {u: k for k, u in enumerate(self.SEQ)}
        self._users = None

    def forward(self, seed):
        v = [0] * NV
        for k, val in seed.items():
            v[k] = val
        ns = {'v': v, '__builtins__': {}}
        for u, i, kd in self.SOLVE:
            _solvevar(v, ns, u, i, kd)
        return v

    def badatoms(self, v):
        ns = {'v': v, '__builtins__': {}}
        out = {}
        for i in range(len(atoms)):
            r = eval(acodes[i], ns)
            if r:
                out[i] = r
        return out

    def seed_of(self, v):
        return {f: v[f] for f in self.FREE if v[f] != 0}

    # ---- incremental ----
    def users(self):
        if self._users is None:
            us = collections.defaultdict(list)
            for w in self.SEQ:
                i, _ = self.definer[w]
                for u in avars[i]:
                    if u != w:
                        us[u].append(w)
            self._users = dict(us)
        return self._users

    def downstream(self, changed):
        us = self.users()
        aff = set(); st = list(changed)
        while st:
            u = st.pop()
            for w in us.get(u, ()):
                if w not in aff:
                    aff.add(w); st.append(w)
        return aff

    def apply_delta(self, v0, changes):
        v = list(v0)
        for k, val in changes.items():
            v[k] = val
        aff = self.downstream(changes.keys())
        ns = {'v': v, '__builtins__': {}}
        for u in sorted(aff, key=lambda u: self._pos[u]):
            i, kind = self.definer[u]
            _solvevar(v, ns, u, i, kind[0])
        return v, aff


_ATOM_OF = collections.defaultdict(list)
for _i, _vs in enumerate(avars):
    for _u in _vs:
        _ATOM_OF[_u].append(_i)


def validate(eng, vd):
    """Demotion must be score-neutral: seeding with current values reproduces vd exactly."""
    v = eng.forward(eng.seed_of(vd))
    return sum(1 for u in range(NV) if v[u] != vd[u])
