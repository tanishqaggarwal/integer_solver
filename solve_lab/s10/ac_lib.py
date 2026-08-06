"""S12 shared fast infrastructure for activation experiments.

* localized forward evaluation (only the downstream cone of the changed seeds)
* incremental cost measurement (touched atoms -> touched equations only)
* cluster gradient support
"""
import os, sys, collections, time
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad
P = ad.P
definer, atom_out = L.definer, L.atom_out
FREE = set(ad.FREE); FORBID = {2081, 4287}
ORDER = ad.ORDER
OPOS = {t: i for i, t in enumerate(ORDER)}
CHECKSET = set(a for a in range(L.NA) if a not in atom_out)

# forward successor graph over defined variables
succ = collections.defaultdict(list)
for _t in definer:
    _a = definer[_t]
    for _w in L.avars[_a]:
        if _w != _t: succ[_w].append(_t)

_dcache = {}
def down(seed):
    """structural downstream cone of a variable (defined vars it can influence)."""
    if seed in _dcache: return _dcache[seed]
    seen, st = set(), [seed]
    while st:
        x = st.pop()
        for t in succ.get(x, ()):
            if t not in seen:
                seen.add(t); st.append(t)
    _dcache[seed] = seen
    return seen

def fwd_local(v, seeds, rounds=6):
    """ad.fwd restricted to the downstream cone of `seeds` -- identical result,
    far cheaper."""
    reg = set()
    for s in seeds: reg |= down(s)
    ords = sorted((OPOS[t] for t in reg if t in OPOS))
    seq = [ORDER[i] for i in ords]
    for _ in range(rounds):
        ch = False
        for u in seq:
            nv = T.solve_lin(definer[u], u, v)
            if nv is not None and nv != v[u]:
                v[u] = nv; ch = True
        if not ch: break
    return v

class Base:
    def __init__(self, path):
        self.v0 = L.load(path)
        self.av0 = L.all_atom_values(self.v0)
        self.nz0 = set(a for a in range(L.NA) if self.av0[a])
        self.fail0 = set(L.failing_eqs(self.av0))
        self.score0 = L.NEQ - len(self.fail0)
        self.vm0 = [x % P for x in self.v0]

    def cost(self, v, changed):
        """(score, newly-nonzero atoms, of which checks, equations newly failing)."""
        ta = L.touched_atoms(v, self.av0, changed)
        nz = set(self.nz0)
        for a, val in ta.items():
            if val: nz.add(a)
            else: nz.discard(a)
        av = self.av0
        if ta:
            av = list(self.av0)
            for a, val in ta.items(): av[a] = val
        eqs = L.eqs_of_atoms(nz)
        fail = set(i for i in eqs if L.eq_value(i, av) != 0)
        newnz = nz - self.nz0
        return (L.NEQ - len(fail), newnz, newnz & CHECKSET,
                fail - self.fail0, self.fail0 - fail, av, nz)

def grad_supp(v, bad):
    vm = [x % P for x in v]
    s = set()
    for a in bad: s |= set(ad.grad(a, vm))
    return s
