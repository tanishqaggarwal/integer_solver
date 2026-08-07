#!/usr/bin/env python3
"""Agent Y -- run the bound machinery at ANY centre D, with no table rebuild.

THE UNIFICATION.  For any centre D subset of {0..255}, with k_D = sum_{i in D} 2^i,

    k  =  k_D  +  sum_{i in S\\D} 2^i  -  sum_{i in D\\S} 2^i

so the residual target  T - fold(D)  is a SIGNED low-weight combination of ladder points
whenever S is close to D.  Exhausting weight M there proves  hamming_distance(S, D) > M.
D = {} is the forward weight bound (agent X); D = {0..255} is the complement bound (agent Y).

WHAT THE UNSIGNED ENGINE CAN AND CANNOT DO -- read this before quoting a result.
The engine in this directory is UNSIGNED: the table holds sums of +2^i G and the scan
subtracts them from a base.  It therefore certifies only the ONE-SIDED balls:

  "up"   target  T - fold(D)   covers  S superset of D  with |S \\ D| <= 9
  "down" target  fold(D) - T   covers  S subset   of D  with |D \\ S| <= 9

The full two-sided Hamming ball around a general D needs a SIGNED table (+-L_i on both
sides) and is NOT what this engine computes.  Note the two centres already swept are exactly
the two where the signed problem degenerates to the unsigned one -- every S contains {} and
every S is contained in {0..255} -- which is a third, independent reason those two centres
are special, alongside the XOR-affinity and affine-self-map uniqueness proofs.

usage:
    python3 ycentre.py 3,17,88,201            # explicit centre
    python3 ycentre.py -                      # empty centre (reproduces agent X's target T)
    python3 ycentre.py all                    # full centre (reproduces T')
then, with no table rebuild:
    ./ymitm scan data_centre_<tag>_up.txt   B ../agentX_work/tbl4s.bin ../agentX_work/bm4.bin rep.txt
    ./ymitm scan data_centre_<tag>_down.txt B ../agentX_work/tbl4s.bin ../agentX_work/bm4.bin rep.txt
for B = 2,3,4,5 (covers distance 3..9; distance <= 4 is the direct table probe, see yedge.py).
"""
import json, os, sys, hashlib

HERE = os.path.dirname(os.path.abspath(__file__))
d = json.load(open(os.path.join(HERE, 'ydata.json')))
p = int(d['p']); A_ = int(d['a']); N = int(d['N'])
L = [(int(x), int(y)) for x, y in d['ladder']]
T = (int(d['T'][0]), int(d['T'][1]))
A = (int(d['A'][0]), int(d['A'][1]))
Tp = (int(d['Tp'][0]), int(d['Tp'][1]))

def inv(z): return pow(z, p - 2, p)
def add(P, Q):
    if P is None: return Q
    if Q is None: return P
    x1, y1 = P; x2, y2 = Q
    if (x1 - x2) % p == 0:
        if (y1 + y2) % p == 0: return None
        l = (3 * x1 * x1 + A_) % p * inv(2 * y1 % p) % p
    else:
        l = (y2 - y1) % p * inv((x2 - x1) % p) % p
    x3 = (l * l - x1 - x2) % p
    return (x3, (l * (x1 - x3) - y1) % p)
def neg(P): return None if P is None else (P[0], (-P[1]) % p)
def sub(P, Q): return add(P, neg(Q))

arg = sys.argv[1]
if arg == '-':      D = []
elif arg == 'all':  D = list(range(256))
else:               D = sorted(set(int(v) for v in arg.split(',')))
assert all(0 <= i < 256 for i in D)

fold_D = None
for i in D: fold_D = add(fold_D, L[i])
k_D = sum(1 << i for i in D)

up   = sub(T, fold_D)          # covers S superset of D
down = sub(fold_D, T)          # covers S subset of D
assert down == neg(up) or (up is None and down is None)

tag = ('empty' if not D else 'all' if len(D) == 256 else
       hashlib.sha256((','.join(map(str, D))).encode()).hexdigest()[:10])
lad = ''.join('%d %d\n' % (x, y) for x, y in L)
made = []
for nm, P in (('up', up), ('down', down)):
    if P is None:
        print('  %-4s target is the identity O -- S == D exactly; nothing to scan' % nm)
        continue
    path = os.path.join(HERE, 'data_centre_%s_%s.txt' % (tag, nm))
    with open(path, 'w') as f:
        f.write('%d %d\n' % (P[0], P[1])); f.write(lad)
    made.append(path)

print('centre D          : |D| = %d  %s' % (len(D), (D if len(D) <= 12 else '(%d indices)' % len(D))))
print('k_D               : %d' % k_D)
print('tag               : %s' % tag)
print('up   = T - fold(D): %s' % (up,))
print('down = fold(D) - T: %s' % (down,))
# self-checks against the two centres already swept
if not D:
    print('CHECK  up == T  (agent X\'s target)            : %s' % (up == T))
if len(D) == 256:
    print('CHECK  fold(D) == A = (2^256-1)G              : %s' % (fold_D == A))
    print('CHECK  down == T\' (this thread\'s target)      : %s' % (down == Tp))
    print('CHECK  k_D == 2^256 - 1                       : %s' % (k_D == 2**256 - 1))
for m in made: print('wrote %s' % m)
