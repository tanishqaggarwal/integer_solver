#!/usr/bin/env python3
"""newton.py -- Newton polytope of the modmul polynomial and of a comb window.

Modmul identity (schoolbook): P = sum_{i,j} 2^{i+j} a_i b_j - sum_k 2^k c_k
                                   - sum_l 2^l q_l * (p-terms) + const.
As a polynomial in the bit variables, the monomials are:
  a_i b_j   (bilinear, the products)
  c_k, q_l  (linear)
Newton polytope = conv hull of exponent vectors. We analyze:
  * Newt(A*B) = Minkowski sum of two simplices (bilinear signature)
  * whether it Minkowski-decomposes (Karatsuba/Toom structure)
  * whether any of this shrinks the {0,1} search space.
"""
import os
import sys
import json
import itertools
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


def newt_AB(s):
    """exponent vectors of A*B monomials a_i b_j in R^{2s} (a-block, b-block)."""
    pts = []
    for i in range(s):
        for j in range(s):
            e = np.zeros(2 * s, dtype=int)
            e[i] = 1
            e[s + j] = 1
            pts.append(e)
    return np.array(pts)


def analyze_hull(pts, label):
    from scipy.spatial import ConvexHull
    # embed into affine hull: subtract mean, compute rank
    P = pts.astype(float)
    c = P.mean(axis=0)
    Q = P - c
    # rank = affine dimension
    u, sv, vt = np.linalg.svd(Q)
    dim = int((sv > 1e-9).sum())
    info = dict(label=label, n_monomials=len(pts), ambient=pts.shape[1],
                affine_dim=dim)
    if dim >= 1:
        # project onto principal axes to a full-dim space for ConvexHull
        basis = vt[:dim].T
        proj = Q @ basis
        if dim == 1:
            info['n_vertices'] = 2
        else:
            try:
                h = ConvexHull(proj)
                info['n_vertices'] = len(np.unique(h.vertices))
                info['n_facets'] = len(h.simplices)
                info['volume_normalized'] = round(h.volume, 6)
            except Exception as e:
                info['hull_error'] = str(e)
    return info


def is_product_of_simplices(s):
    """Newt(A*B) should be combinatorially Delta_{s-1} x Delta_{s-1}:
       affine dim 2(s-1), (s)^2 vertices? No -- product of two (s-1)-simplices
       has s*s vertices and dimension 2(s-1). Check vertex count and dim."""
    pts = newt_AB(s)
    info = analyze_hull(pts, f"Newt(A*B) s={s}")
    expected_dim = 2 * (s - 1)
    # every a_i b_j is a vertex of Delta x Delta
    info['expected_dim_Delta_x_Delta'] = expected_dim
    info['expected_vertices'] = s * s
    info['matches_product_of_simplices'] = (
        info['affine_dim'] == expected_dim and
        info.get('n_vertices', -1) == s * s)
    return info


def minkowski_check(s):
    """Newt(A*B) = Newt(A) (+) Newt(B): the Minkowski sum of two standard
    simplices. Verify the vertex set equals the Minkowski sum of the two
    simplices' vertices."""
    simplexA = [tuple(1 if k == i else 0 for k in range(s)) for i in range(s)]
    simplexB = [tuple(1 if k == j else 0 for k in range(s)) for j in range(s)]
    mink = set()
    for va in simplexA:
        for vb in simplexB:
            mink.add(va + vb)
    actual = set(tuple(e) for e in newt_AB(s))
    return dict(minkowski_sum_size=len(mink), monomial_set_size=len(actual),
                equal=(mink == actual))


def comb_window_atoms(w):
    """A comb window: one-hot selector u_0..u_{D-1}, D=2^w; the selected atom is
    sum_t u_t * table[t]. As a polynomial the monomials are the D linear terms
    u_t. Newton polytope = standard simplex Delta_{D-1}."""
    D = 1 << w
    pts = np.eye(D, dtype=int)
    info = analyze_hull(pts, f"comb window w={w} (D={D})")
    info['is_standard_simplex'] = (info['affine_dim'] == D - 1)
    return info


if __name__ == '__main__':
    out = {}
    out['product_of_simplices'] = [is_product_of_simplices(s) for s in [2, 3, 4, 5]]
    out['minkowski'] = [dict(s=s, **minkowski_check(s)) for s in [2, 3, 4, 5]]
    out['comb_windows'] = [comb_window_atoms(w) for w in [1, 2, 3, 4]]
    print(json.dumps(out, indent=1, default=str))
    with open('newton.json', 'w') as f:
        json.dump(out, f, indent=1, default=str)
