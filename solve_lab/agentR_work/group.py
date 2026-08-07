#!/usr/bin/env python3
"""Exact group arithmetic for the fold law, in shifted coordinates Y^2 = X^3 + B (measured: a=0)."""
from model import P, K, S, inv, to_short, from_short, load_points, TARGET

def fit_B(pts_short):
    X, Y = pts_short[0]
    return (Y * Y - X ** 3) % P

def add(A, Bp):
    if A is None: return Bp
    if Bp is None: return A
    ax, ay = A; bx, by = Bp
    if ax == bx:
        if (ay + by) % P == 0: return None
        l = 3 * ax * ax % P * inv(2 * ay) % P
    else:
        l = (by - ay) % P * inv(bx - ax) % P
    cx = (l * l - ax - bx) % P
    return (cx, (l * (ax - cx) - ay) % P)

def neg(A):
    return None if A is None else (A[0], (-A[1]) % P)

def mul(k, A):
    R = None; Q = A
    if k < 0: k, Q = -k, neg(A)
    while k:
        if k & 1: R = add(R, Q)
        Q = add(Q, Q); k >>= 1
    return R
