#!/usr/bin/env python3
"""ACCUMULATOR TEST.

With two same-tree booleans ON, what value does the selected wire pair take?
Hypotheses tested, by freezing (x1,y1) at the candidate value and running the chain repair:
  sum   : coordinate-wise sum mod p of the two booleans' forced pairs
  chord : the instance's own degree-3 combination law  x3 = l^2 - x1 - x2 - K , y3 = l(x1-x3) - y1
  first / second : one boolean's pair alone
  zero  : what the greedy repair picks on its own
Reports the number of nonzero residual atoms and the checker-equivalent score for each.
"""
import sys, os, json, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import gs2
from fwd import NV
E = gs2.E
p = 115792089237316195423570985008687907853269984665640564039457584007908834671663
K = 97553848499418123410591666447050222001188385549510401465815187079080512838891
K1 = 91416258160755509149180373473728639746431157665678710450404458852172057265575180278101002
K2 = 125787314747601108116039725163361763116550465675981151838811516827327919228823597744635626
X1, Y1, X2, Y2, X3, Y3 = 12186, 16742, 14853, 24908, 22162, 30213


def chordK(P, Q):
    x1, y1 = P; x2, y2 = Q
    if (x2 - x1) % p == 0:
        return None
    l = (y2 - y1) * pow(x2 - x1, p - 2, p) % p
    x3 = (l * l - x1 - x2 - K) % p
    return (x3, (l * (x1 - x3) - y1) % p)


def run(bits, freeze=None, extra=None):
    v = [0] * NV
    for k in bits:
        v[k] = 1
    for k, x in {X3: K1, Y3: K2, 24468: K1, 18956: K2}.items():
        v[k] = x
    fr = set(bits) | {X3, Y3, 24468, 18956}
    if freeze:
        for k, x in freeze.items():
            v[k] = x; fr.add(k)
    if extra:
        for k, x in extra.items():
            v[k] = x
    v, ok = gs2.solve(v, verbose=False, frozen=fr)
    r = E.run(v)
    bad = E.score(r)
    nz = [i for i in range(len(r)) if r[i]]
    return v, 39033 - len(bad), nz


def forced_pairs():
    d = json.load(open(os.path.join(HERE, 'sweep_ii.json')))
    out = {}
    for k, val in d.items():
        f = val.get('forced')
        if f and None not in f:
            out[int(k)] = (int(f[0]), int(f[1]))
    return out


if __name__ == '__main__':
    FP = forced_pairs()
    partner = 5090
    cands = [b for b in FP if b != partner]
    j1, j2 = sorted(cands)[:2]
    P1, P2 = FP[j1], FP[j2]
    print('tree-A booleans %d and %d' % (j1, j2), flush=True)
    print('  P(%d) = %s' % (j1, [str(z)[:30] for z in P1]), flush=True)
    print('  P(%d) = %s' % (j2, [str(z)[:30] for z in P2]), flush=True)
    hyp = {
        'zero':   None,
        'first':  P1,
        'second': P2,
        'sum':    ((P1[0] + P2[0]) % p, (P1[1] + P2[1]) % p),
        'chord':  chordK(P1, P2),
    }
    base, sb, nzb = run([j1, partner])
    print('CONTROL single boolean %d : score %d, nonzero atoms %d' % (j1, sb, len(nzb)), flush=True)
    for name, val in hyp.items():
        t0 = time.time()
        if val is None:
            v, s, nz = run([j1, j2, partner])
        else:
            v, s, nz = run([j1, j2, partner], freeze={X1: val[0], Y1: val[1]})
        print('  %-7s -> score %d, nonzero atoms %d   (x1=%s)  t=%.0fs'
              % (name, s, len(nz), str(v[X1] % p)[:26], time.time() - t0), flush=True)
        for i in nz[:8]:
            print('        NZ', E.res[i][:88], flush=True)
