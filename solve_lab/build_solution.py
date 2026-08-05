#!/usr/bin/env python3
"""Regenerate the complete integer solution from best_partial_39019 in one pass.

The two nonzero atoms at best_partial are product-slacks that fix the twist:
  H: x_9982  = x_12518*x_9897  = -(x_17728 - x_3183)          (absorbs invariant gap G)
  F: x_26977 = x_20510*x_31302 = 6033033*(x_18274 - x_9770)
x_12518 and x_20510 are hubs in the 220-var identity 'wire' (all = +/-x_15); the wire is
EXACTLY quiet (V^1 and V^2+ coefficients vanish in every atom), so it is a free parameter.
Set the wire to sign*V (V=1) and the two rare partners x_9897=-G, x_31302=F0. Done.
"""
import json, sys
from propagate import load_atoms

def main():
    A = load_atoms()
    best = json.load(open('best/best_partial_39019.json'))
    v = {int(k[2:]): val for k, val in best.items()}
    g = lambda i: v.get(i, 0)

    # signed identity union-find: x = sign * rep, over 2-term |c1|==|c2| atoms
    parent, sgn = {}, {}
    def find(x):
        parent.setdefault(x, x); sgn.setdefault(x, 1)
        if parent[x] == x: return x, 1
        r, s = find(parent[x]); parent[x] = r; sgn[x] = sgn[x] * s; return r, sgn[x]
    def union(a, b, rel):
        ra, sa = find(a); rb, sb = find(b)
        if ra != rb: parent[ra] = rb; sgn[ra] = rel * sb * sa
    for poly in A:
        if len(poly) == 2:
            (m1, c1), (m2, c2) = list(poly.items())
            if len(m1) == 1 and len(m2) == 1 and abs(c1) == abs(c2):
                union(m1[0], m2[0], -1 if (c1 > 0) == (c2 > 0) else 1)
    r15 = find(15)[0]
    wire = {y: find(y)[1] for y in list(parent) if find(y)[0] == r15}

    G  = g(17728) - g(3183)
    F0 = 6033033 * (g(18274) - g(9770))
    V = 1
    cand = dict(v)
    for w, s in wire.items():
        cand[w] = s * V
    cand[9982]  = -G
    cand[26977] = F0
    cand[9897]  = -G // (wire[12518] * V)   # x_12518 = sign*V
    cand[31302] = F0 // (wire[20510] * V)   # x_20510 = sign*V

    # sanity on the four fixed atoms
    assert cand[9982]  - cand[9897] * cand[12518] == 0        # atom 1818
    assert cand[26977] - cand[20510] * cand[31302] == 0       # atom 1816
    assert (cand[17728] - cand[3183]) + cand[9982] == 0       # H
    assert 6033033 * (cand[18274] - cand[9770]) - cand[26977] == 0  # F

    out = {f'x_{i}': cand.get(i, 0) for i in range(38748)}
    path = sys.argv[1] if len(sys.argv) > 1 else 'best/SOLUTION.json'
    json.dump(out, open(path, 'w'))
    print(f'wrote {path}: {len(out)} vars, {sum(1 for x in out.values() if x)} nonzero')

if __name__ == '__main__':
    main()
