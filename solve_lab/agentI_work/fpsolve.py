#!/usr/bin/env python3
"""Solve the atom system mod p by propagation + one-at-a-time free-input choice.

Usage: python3 fpsolve.py <policy> [seed]
  policy: wit | zero | rand
"""
import pickle, os, collections, sys, random, time
from model import Model, load_assign
from fp import FpEngine, P

HERE = os.path.dirname(os.path.abspath(__file__))
NV = 38748


def run(policy='wit', seed=1, verbose=True):
    M = Model(); E = FpEngine(M)
    wit = load_assign(os.path.join(HERE, '..', 'best',
                                   'new_instance_partial_39026.json'))
    rng = random.Random(seed)
    val = [None] * NV
    decisions = []
    t0 = time.time()
    conflicts_total = []
    while True:
        n, conf, br = E.propagate(val)
        if conf:
            conflicts_total.extend(conf)
            if verbose:
                print("CONFLICT", collections.Counter(k for k, _ in conf))
                for k, a in conf[:5]:
                    print("   ", k, a, M.src[a][:120])
            break
        if br:
            for u, a, roots in br:
                if val[u] is not None:
                    continue
                if policy == 'wit':
                    w = wit[u] % P
                    val[u] = w if w in roots else roots[0]
                elif policy == 'zero':
                    val[u] = 0 if 0 in roots else roots[0]
                else:
                    val[u] = rng.choice(roots)
                decisions.append((u, 'branch'))
            continue
        # stuck: designate a free input
        unknown = [v for v in range(NV) if val[v] is None]
        if not unknown:
            break
        deg = collections.Counter()
        for a in range(M.na):
            miss = [x for x in E.avarlist[a] if val[x] is None]
            if len(miss) == 2:
                for x in miss:
                    deg[x] += 1
        if deg:
            u = deg.most_common(1)[0][0]
        else:
            u = unknown[0]
        if policy == 'wit':
            val[u] = wit[u] % P
        elif policy == 'zero':
            val[u] = 0
        else:
            val[u] = rng.randrange(P)
        decisions.append((u, 'free'))
        if verbose and len(decisions) % 200 == 0:
            print(f"  decisions={len(decisions)} known={sum(1 for x in val if x is not None)} "
                  f"t={time.time()-t0:.0f}s", flush=True)
    known = sum(1 for x in val if x is not None)
    filled = [0 if x is None else x for x in val]
    bad = [a for a in range(M.na) if E.eval_atom(a, filled) != 0]
    print(f"policy={policy} known={known}/{NV} decisions={len(decisions)} "
          f"conflicts={len(conflicts_total)} atoms_bad_modp={len(bad)} "
          f"t={time.time()-t0:.0f}s")
    return val, filled, bad, conflicts_total, M, E


if __name__ == '__main__':
    pol = sys.argv[1] if len(sys.argv) > 1 else 'wit'
    sd = int(sys.argv[2]) if len(sys.argv) > 2 else 1
    val, filled, bad, conf, M, E = run(pol, sd)
    pickle.dump(val, open(os.path.join(HERE, f'fpsolve_{pol}_{sd}.pkl'), 'wb'))
    for a in bad[:30]:
        print("   bad a%d %s" % (a, M.src[a][:130]))
