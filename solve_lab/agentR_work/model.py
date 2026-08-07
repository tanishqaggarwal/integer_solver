#!/usr/bin/env python3
"""Agent R: reduced-problem model.

Loads, from agent F's decode artifacts (read-only), the reduced instance:
  * prime p, offset K, the degree-3 stage law,
  * the 256 conditional-pin booleans and their forced coordinate pairs,
  * the target root pair.
Provides exact F_p arithmetic for the law and the derived fold operation.
No curve/PRNG framing: everything here is a measured consequence of the law in
stage_roles.json, which every one of F's 72 full stages satisfies.
"""
import os, json, pickle, sys

F = '/home/user/integer_solver/solve_lab/agentF_work'
P = 115792089237316195423570985008687907853269984665640564039457584007908834671663
K = 97553848499418123410591666447050222001188385549510401465815187079080512838891
# root target: the values the instance demands at the root pair (x22162, x30213)
T1 = 91416258160755509149180373473728639746431157665678710450404458852172057265575180278101002
T2 = 125787314747601108116039725163361763116550465675981151838811516827327919228823597744635626
TARGET = (T1 % P, T2 % P)

INF = None   # identity of the fold


def inv(a):
    return pow(a % P, P - 2, P)


def chordK(A, B):
    """The stage law, exactly as stage_law2.py verified it on 72/72 stages."""
    if A is INF: return B
    if B is INF: return A
    ax, ay = A; bx, by = B
    if ax == bx:
        if (ay + by) % P == 0:
            return INF
        l = (3 * ax * ax + 2 * K * ax) % P * inv(2 * ay) % P   # see doc/tangent note
    else:
        l = (by - ay) % P * inv(bx - ax) % P
    ox = (l * l - ax - bx - K) % P
    oy = (l * (ax - ox) - ay) % P
    return (ox, oy)


# ---- shifted coordinates: X = x + K/3 removes the offset entirely ----
S = K * inv(3) % P


def to_short(Pt):
    if Pt is INF: return INF
    return ((Pt[0] + S) % P, Pt[1] % P)


def from_short(Pt):
    if Pt is INF: return INF
    return ((Pt[0] - S) % P, Pt[1] % P)


def load_points():
    """256 booleans -> (x, y) in the instance's coordinates."""
    sw = json.load(open(os.path.join(F, 'sweep_ii.json')))
    pins = json.load(open(os.path.join(F, 'pins.json')))
    pts, tree, missing = {}, {}, []
    for b, r in sw.items():
        tree[int(b)] = r.get('tree')
        if r.get('all_forced'):
            pts[int(b)] = (int(r['forced'][0]) % P, int(r['forced'][1]) % P)
        else:
            missing.append(int(b))
    return pts, tree, missing, pins


if __name__ == '__main__':
    pts, tree, missing, pins = load_points()
    print('points forced:', len(pts), ' missing:', sorted(missing))
    print('tree A:', sum(1 for v in tree.values() if v == 'A'),
          ' tree B:', sum(1 for v in tree.values() if v == 'B'))
    print('target =', TARGET)
