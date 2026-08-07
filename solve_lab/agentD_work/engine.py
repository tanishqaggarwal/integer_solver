"""Agent D search engine.

Invariant of the search space: all 31,475 gate atoms are zero (we always forward-
evaluate).  State = assignment of the 7,273 free inputs.  Objective = number of
satisfied equations, computed exactly over Z.

Core operations:
  * forward(v)                  full topological gate evaluation
  * cone_free(target_vars)      free inputs that can influence the targets
  * probe_linear(v, u, c)       exact (c0, slope, is_linear) of check atom c in free input u
  * repair(v, c)                candidate moves closing check c
"""
import collections, heapq, random, time, sys
import dlib as L

P = L.P
NV = L.NVARS

# ---------- influence graph (var -> vars directly computed from it) ----------
_infl = collections.defaultdict(list)
for _t, _a in L.definer.items():
    for _u in L.avars[_a]:
        if _u != _t:
            _infl[_u].append(_t)
infl = dict(_infl)


def cone_free(targets, maxn=None):
    """Free inputs that can reach any of `targets` (variables) through gate edges."""
    # reverse: build via backward BFS using definer
    seen = set()
    stack = list(targets)
    out = set()
    while stack:
        u = stack.pop()
        if u in seen:
            continue
        seen.add(u)
        a = L.definer.get(u)
        if a is None:
            out.add(u)
            continue
        for w in L.avars[a]:
            if w != u and w not in seen:
                stack.append(w)
    return out, seen


def forward(v, block=()):
    block = set(block)
    for t in L.topo:
        a = L.definer[t]
        if a in block:
            continue
        nv = L.solve_for(a, t, v)
        if nv is not None:
            v[t] = nv


def fwd_from(v, seeds, block=()):
    """Incremental forward: only recompute vars downstream of `seeds`."""
    block = set(block)
    dirty = set()
    pq = []
    for u in seeds:
        for w in infl.get(u, ()):
            if w not in dirty:
                dirty.add(w)
                heapq.heappush(pq, (L.tidx.get(w, L.BIG), w))
    changed = set(seeds)
    while pq:
        _, t = heapq.heappop(pq)
        dirty.discard(t)
        a = L.definer[t]
        if a in block:
            continue
        nv = L.solve_for(a, t, v)
        if nv is None or nv == v[t]:
            continue
        v[t] = nv
        changed.add(t)
        for w in infl.get(t, ()):
            if w not in dirty:
                dirty.add(w)
                heapq.heappush(pq, (L.tidx.get(w, L.BIG), w))
    return changed


def score_full(v):
    av = L.all_atom_values(v)
    f = L.failing_eqs(av)
    return L.NEQ - len(f), av, f


def nz_checks(av):
    return [a for a in range(L.NA) if av[a]]


# ---------- scoring cache ----------
class State:
    def __init__(self, v):
        self.v = list(v)
        self.refresh()

    def refresh(self):
        self.av = L.all_atom_values(self.v)
        self.fail = set(L.failing_eqs(self.av))
        self.score = L.NEQ - len(self.fail)

    def copy(self):
        s = State.__new__(State)
        s.v = list(self.v)
        s.av = list(self.av)
        s.fail = set(self.fail)
        s.score = self.score
        return s

    def apply(self, seeds, block=()):
        """Set free inputs, incrementally forward, incrementally rescore."""
        for u, val in seeds.items():
            self.v[u] = val
        changed = fwd_from(self.v, list(seeds), block=block)
        cand = set()
        for u in changed:
            cand.update(L.var_atoms.get(u, ()))
        eqs = set()
        for a in cand:
            nv = L.evalpoly(L.polys[a], self.v)
            if nv != self.av[a]:
                self.av[a] = nv
                eqs.update(L.atom2eq.get(a, ()))
        for i in eqs:
            m, sq, co = L.eq_atoms[i]
            s = 0
            for a, c in co.items():
                s += c * self.av[a]
            if s:
                self.fail.add(i)
            else:
                self.fail.discard(i)
        self.score = L.NEQ - len(self.fail)
        return changed
