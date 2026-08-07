#!/usr/bin/env python3
"""Lift a mod-p solution to Z: seed the free inputs with their residues and
propagate exactly over Z.  Reports where the lift breaks."""
import pickle, os, sys, collections, json, time
from model import Model, load_assign
from prop import Engine
HERE = os.path.dirname(os.path.abspath(__file__))
P = 2**256 - 2**32 - 977
NV = 38748


def main(tag):
    d = pickle.load(open(os.path.join(HERE, f'fprun_{tag}.pkl'), 'rb'))
    valp = d['val']; decisions = set(d['decisions'])
    M = Model(); E = Engine(M)
    z = [None] * NV
    for v in decisions:
        if valp[v] is not None:
            z[v] = valp[v]
    n, conf, br = E.propagate(z)
    known = sum(1 for x in z if x is not None)
    print(f"lift: seeded {len(decisions)} decisions -> known {known}/{NV}, "
          f"conflicts {len(conf)}, branch {len(br)}")
    cc = collections.Counter(k for k, _ in conf)
    print("conflict kinds:", cc)
    for k, a in conf[:20]:
        print("   ", k, a, M.src[a][:130])
    # consistency with mod p
    dis = [v for v in range(NV) if z[v] is not None and valp[v] is not None
           and z[v] % P != valp[v]]
    print("vars whose Z lift disagrees with mod-p value:", len(dis), dis[:20])
    zz = [0 if x is None else x for x in z]
    fails, av, cv = M.eq_fail(zz)
    print("SCORE of naive lift:", M.ne - len(fails), "/", M.ne)
    json.dump({f"x_{i}": zz[i] for i in range(NV) if zz[i]},
              open(os.path.join(HERE, f'lift_{tag}.json'), 'w'))
    return z


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else 'wit_1')
