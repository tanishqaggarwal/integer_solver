"""State with undo + repair move generation."""
import collections, heapq, random, time, sys
import dlib as L
import rad

P = L.P

_infl = collections.defaultdict(list)
for _t, _a in L.definer.items():
    for _u in L.avars[_a]:
        if _u != _t:
            _infl[_u].append(_t)
infl = dict(_infl)


class St:
    def __init__(self, v):
        self.v = list(v)
        self.av = L.all_atom_values(self.v)
        self.fail = set(L.failing_eqs(self.av))
        self.score = L.NEQ - len(self.fail)

    def clone(self):
        s = St.__new__(St)
        s.v = list(self.v)
        s.av = list(self.av)
        s.fail = set(self.fail)
        s.score = self.score
        return s

    def apply(self, seeds):
        """Apply free-input changes; returns undo record."""
        oldv = {}
        pq = []
        dirty = set()
        for u, val in seeds.items():
            if self.v[u] == val:
                continue
            oldv[u] = self.v[u]
            self.v[u] = val
            for w in infl.get(u, ()):
                if w not in dirty:
                    dirty.add(w)
                    heapq.heappush(pq, (L.tidx.get(w, L.BIG), w))
        while pq:
            _, t = heapq.heappop(pq)
            dirty.discard(t)
            a = L.definer[t]
            nv = L.solve_for(a, t, self.v)
            if nv is None or nv == self.v[t]:
                continue
            if t not in oldv:
                oldv[t] = self.v[t]
            self.v[t] = nv
            for w in infl.get(t, ()):
                if w not in dirty:
                    dirty.add(w)
                    heapq.heappush(pq, (L.tidx.get(w, L.BIG), w))
        cand = set()
        for u in oldv:
            cand.update(L.var_atoms.get(u, ()))
        olda = {}
        eqs = set()
        for a in cand:
            nv = L.evalpoly(L.polys[a], self.v)
            if nv != self.av[a]:
                olda[a] = self.av[a]
                self.av[a] = nv
                eqs.update(L.atom2eq.get(a, ()))
        addf = []
        delf = []
        for i in eqs:
            m, sq, co = L.eq_atoms[i]
            s = 0
            for a, c in co.items():
                s += c * self.av[a]
            if s:
                if i not in self.fail:
                    self.fail.add(i)
                    addf.append(i)
            else:
                if i in self.fail:
                    self.fail.discard(i)
                    delf.append(i)
        self.score = L.NEQ - len(self.fail)
        return (oldv, olda, addf, delf)

    def revert(self, rec):
        oldv, olda, addf, delf = rec
        for u, val in oldv.items():
            self.v[u] = val
        for a, val in olda.items():
            self.av[a] = val
        for i in addf:
            self.fail.discard(i)
        for i in delf:
            self.fail.add(i)
        self.score = L.NEQ - len(self.fail)

    def nz(self):
        return [a for a in range(L.NA) if self.av[a]]


def forward_all(v):
    for t in L.topo:
        a = L.definer[t]
        nv = L.solve_for(a, t, v)
        if nv is not None:
            v[t] = nv


def probe_solve(st, c, u):
    """Exact: is atom c affine in free input u?  If so return the required new value
    of u making c zero, else None."""
    c0 = st.av[c]
    if c0 == 0:
        return None
    base = st.v[u]
    r1 = st.apply({u: base + 1})
    c1 = st.av[c]
    st.revert(r1)
    if c1 == c0:
        return None
    r2 = st.apply({u: base + 2})
    c2 = st.av[c]
    st.revert(r2)
    if c2 - c1 != c1 - c0:
        return None                      # not affine along u
    slope = c1 - c0
    if c0 % slope:
        return None
    return base - c0 // slope


def repair_candidates(st, c, knobs=None, limit=None, shuffle=False):
    """Return list of (u, newval) that exactly zero atom c."""
    if knobs is None:
        knobs = sorted(rad.free_knobs(c, st.v))
    ks = list(knobs)
    if shuffle:
        random.shuffle(ks)
    if limit:
        ks = ks[:limit]
    out = []
    for u in ks:
        nv = probe_solve(st, c, u)
        if nv is not None and nv != st.v[u]:
            out.append((u, nv))
    return out
