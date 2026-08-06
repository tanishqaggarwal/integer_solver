"""bl_core: shared frame machinery + cheap pre-filter + enriched engine wrapper."""
import os, sys, json, time, collections
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad
P = 2**256 - 2**32 - 977
NA, NEQ = L.NA, L.NEQ

# ---------------- booleans ----------------
BOOLATOM = {}
for _a, _p in enumerate(L.polys):
    ks = list(_p.items())
    if len(ks) == 2:
        sq = [m for m, c in ks if len(m) == 2 and m[0] == m[1]]
        li = [m for m, c in ks if len(m) == 1]
        if sq and li and sq[0][0] == li[0][0]: BOOLATOM[li[0][0]] = _a
BOOL = set(BOOLATOM)

FORBID = {2081, 4287}


class Frame:
    """A choice of which gate atoms stay gates.  DETACH maps var->atom to release."""
    def __init__(self, detach=()):
        self.DETACH = dict(detach)
        self.definer = {t: a for t, a in L.definer.items() if t not in self.DETACH}
        self.atom_out = {a: o for a, o in L.atom_out.items()
                         if a not in set(self.DETACH.values())}
        self.ORDER = [t for t in ad.ORDER if t not in self.DETACH]
        self.FREE = set(t for t in range(L.NVARS) if t not in self.definer)
        self.CHECKS = [a for a in range(NA) if a not in self.atom_out]

    def fwd(self, v, rounds=6):
        de, so = self.definer, T.solve_lin
        for _ in range(rounds):
            for u in self.ORDER:
                nv = so(de[u], u, v)
                if nv is not None: v[u] = nv
        return v

    def cone(self, atoms):
        c, st = set(), []
        for a in atoms: st += list(L.avars[a])
        while st:
            t = st.pop()
            if t in c: continue
            c.add(t)
            d = self.definer.get(t)
            if d is None: continue
            for w in L.avars[d]:
                if w != t: st.append(w)
        return c


CANON = Frame()
F2 = Frame({7068: 22229, 28730: 22230, 29854: 35758, 31864: 35761, 642: 35762})


def pot(v):
    av = L.all_atom_values(v)
    nz = [a for a in range(NA) if av[a]]
    s = NEQ - len(L.failing_eqs(av))
    return (s, -len(nz), -sum(abs(av[a]).bit_length() for a in nz)), av, nz


def cheap(v0, flips, F, rounds=6):
    """flip the given vars, forward-solve in frame F, return ((score,-nz,-bits), v)."""
    v = list(v0)
    for u in flips: v[u] = 1 - v[u]
    F.fwd(v, rounds=rounds)
    p, av, nz = pot(v)
    return p, v, nz


# ---------------- engine, parameterised by frame ----------------
def moves(a, v, av, F, nnewton=20):
    out = []
    for w in sorted(set(L.avars[a])):
        if w in FORBID: continue
        tgt = T.solve_lin(a, w, v)
        if tgt is None or tgt == v[w]: continue
        if w in F.FREE: out.append((w, tgt))
        else:
            d = F.definer.get(w)
            if d is None: continue
            vv = list(v); vv[w] = tgt
            for u in sorted(set(L.avars[d])):
                if u == w or u not in F.FREE or u in FORBID: continue
                nv = T.solve_lin(d, u, vv)
                if nv is not None: out.append((u, nv))
    r = av[a] % P
    if r:
        vm = [x % P for x in v]
        try: g = ad.grad(a, vm)
        except Exception: g = {}
        cand = sorted((len(L.var_atoms[u]), u, d) for u, d in g.items()
                      if u not in FORBID and d % P)
        for _, u, d in cand[:nnewton]:
            out.append((u, v[u] + (-r * pow(d, -1, P)) % P))
    return out


def engine(v, F, iters=60, budget=240, nnewton=20, verbose=False, tag=''):
    cur, av, nz = pot(v)
    t0 = time.time()
    for it in range(iters):
        if time.time() - t0 > budget: break
        got = None
        for a in nz:
            for u, nv in moves(a, v, av, F, nnewton=nnewton):
                tr = list(v); tr[u] = nv
                F.fwd(tr, rounds=6)
                p2, av2, nz2 = pot(tr)
                if p2 > cur: got = (a, u, p2, tr, av2, nz2); break
            if got: break
        if not got: break
        a, u, p2, tr, av2, nz2 = got
        if verbose:
            print(f'   {tag} it{it}: a{a} via x_{u}  {cur[0]} -> {p2[0]}  nz {len(nz)}->{len(nz2)}', flush=True)
        v, cur, av, nz = tr, p2, av2, nz2
    return cur, v, nz
