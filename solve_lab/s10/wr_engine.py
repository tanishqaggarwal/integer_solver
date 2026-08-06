"""WR: the enriched repair engine, frame-aware (works inside wr_frame.F_WIRE)."""
import os, sys, time, collections
HERE = '/home/user/integer_solver/solve_lab/s10'
LAB = '/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
import wr_frame as W
P = ad.P


class Engine:
    def __init__(self, F, forbid=()):
        self.F = F
        self.FORBID = set(forbid)

    def pot(self, v):
        av = L.all_atom_values(v)
        nz = [a for a in range(L.NA) if av[a]]
        s = L.NEQ - len(L.failing_eqs(av))
        return (s, -len(nz), -sum(abs(av[a]).bit_length() for a in nz)), av, nz

    def grad(self, check, vm):
        F = self.F
        lam = collections.defaultdict(int)
        for x in L.avars[check]:
            lam[x] = (lam[x] + ad.dpart(check, x, vm)) % P
        for t in reversed(F.ORDER):
            lt = lam.get(t, 0)
            if not lt:
                continue
            a = F.definer[t]
            d = ad.dpart(a, t, vm)
            if d % P == 0:
                continue
            f = -lt * pow(d, -1, P) % P
            for x in L.avars[a]:
                if x == t:
                    continue
                dx = ad.dpart(a, x, vm)
                if dx:
                    lam[x] = (lam[x] + f * dx) % P
            lam[t] = 0
        return {u: lam[u] % P for u in F.FREE if lam.get(u, 0) % P}

    def moves(self, a, v, av, nnewton=25):
        F = self.F
        out = []
        for w in sorted(set(L.avars[a])):
            if w in self.FORBID:
                continue
            tgt = T.solve_lin(a, w, v)
            if tgt is None or tgt == v[w]:
                continue
            if w in F.FREE:
                out.append((w, tgt))
            else:
                d = F.definer.get(w)
                if d is None:
                    continue
                vv = list(v); vv[w] = tgt
                for u in sorted(set(L.avars[d])):
                    if u == w or u not in F.FREE or u in self.FORBID:
                        continue
                    nv = T.solve_lin(d, u, vv)
                    if nv is not None:
                        out.append((u, nv))
        r = av[a] % P
        if r:
            vm = [x % P for x in v]
            try:
                g = self.grad(a, vm)
            except Exception:
                g = {}
            cand = sorted((len(L.var_atoms[u]), u, d) for u, d in g.items()
                          if u not in self.FORBID and d % P)
            for _, u, d in cand[:nnewton]:
                out.append((u, v[u] + (-r * pow(d, -1, P)) % P))
        return out

    def run(self, v, tag, iters=80, budget=2400, save=True):
        F = self.F
        cur, av, nz = self.pot(v)
        print(f'{tag}: start {cur[0]} (nonzero {len(nz)})', flush=True)
        t0 = time.time()
        for it in range(iters):
            if time.time() - t0 > budget:
                print(f'  {tag}: budget out at {cur[0]}', flush=True)
                break
            got = None
            for a in nz:
                for u, nv in self.moves(a, v, av):
                    tr = list(v); tr[u] = nv
                    F.fwd(tr, rounds=6)
                    p2, av2, nz2 = self.pot(tr)
                    if p2 > cur:
                        got = (a, u, p2, tr, av2, nz2)
                        break
                if got:
                    break
            if not got:
                print(f'  {tag} it{it}: STUCK score {cur[0]} nonzero {sorted(nz)}', flush=True)
                break
            a, u, p2, tr, av2, nz2 = got
            print(f'  {tag} it{it}: a{a} via x_{u}  {cur[0]} -> {p2[0]}  '
                  f'nz {len(nz)} -> {len(nz2)}', flush=True)
            v, cur, av, nz = tr, p2, av2, nz2
            if save and p2[0] > 39026:
                T.save(v, os.path.join(HERE, 'wr_best.json'))
                print('   *** saved wr_best.json', flush=True)
        if save:
            T.save(v, os.path.join(HERE, f'wr_engine_{tag}_{cur[0]}.json'))
        print(f'{tag} FINAL {cur[0]} nonzero {sorted(nz)}', flush=True)
        return v, cur


if __name__ == '__main__':
    path = sys.argv[1]
    v = L.load(path if os.path.isabs(path) else os.path.join(HERE, path))
    E = Engine(W.F_WIRE)
    E.run(v, os.path.basename(path).replace('.json', ''),
          budget=int(sys.argv[2]) if len(sys.argv) > 2 else 2400)
