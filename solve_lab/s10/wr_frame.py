"""WR shared frame: frame3's detachments PLUS the wire root x_26064.

In this frame a37694 (the bare wire pin) is a CHECK, so the wire may carry any
value w; the rest of the frame is exactly frame3's, in which the delivered
39,026 witness is on-manifold.
"""
import os, sys, collections
HERE = '/home/user/integer_solver/solve_lab/s10'
LAB = '/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
P = ad.P
ROOT, ROOTATOM = 26064, 37694

F3 = {7068: 22229, 28730: 22230, 29854: 35758, 31864: 35761, 642: 35762, 4432: 22231}


def build(detach):
    definer = {t: a for t, a in L.definer.items() if t not in detach}
    atom_out = {a: o for a, o in L.atom_out.items() if a not in set(detach.values())}
    ORDER = [t for t in ad.ORDER if t not in detach]
    FREE = set(t for t in range(L.NVARS) if t not in definer)
    CHECKS = [a for a in range(L.NA) if a not in atom_out]
    return definer, atom_out, ORDER, FREE, CHECKS


class Frame:
    def __init__(self, detach):
        self.detach = dict(detach)
        (self.definer, self.atom_out, self.ORDER,
         self.FREE, self.CHECKS) = build(self.detach)

    def fwd(self, v, rounds=8):
        for _ in range(rounds):
            for u in self.ORDER:
                nv = T.solve_lin(self.definer[u], u, v)
                if nv is not None:
                    v[u] = nv
        return v

    def report(self, v, tag='', quiet=False):
        av = L.all_atom_values(v)
        nz = [a for a in range(L.NA) if av[a]]
        fail = L.failing_eqs(av)
        sc = L.NEQ - len(fail)
        if not quiet:
            print(f'{tag}: score={sc} failing={len(fail)} nonzero={len(nz)} {sorted(nz)}',
                  flush=True)
        return av, nz, fail, sc


F_CANON = Frame({})
F_F3 = Frame(F3)
F_WIRE = Frame({**F3, ROOT: ROOTATOM})
F_ROOTONLY = Frame({ROOT: ROOTATOM})


def wire_of(base):
    return sorted(u for u in range(L.NVARS) if base[u] == P)


if __name__ == '__main__':
    base = L.load(os.path.join(LAB, 'best', 'new_instance_partial_39026.json'))
    for nm, F in (('canonical', F_CANON), ('frame3', F_F3),
                  ('frame3+root', F_WIRE), ('root only', F_ROOTONLY)):
        v = list(base); F.fwd(v)
        F.report(v, f'deliverable in {nm}')
