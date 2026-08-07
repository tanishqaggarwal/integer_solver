#!/usr/bin/env python3
"""Batched affine group arithmetic (Montgomery's trick) for the fold law's group."""
from model import P

def batch_inv(xs):
    n = len(xs); pre = [1] * (n + 1)
    for i, x in enumerate(xs): pre[i + 1] = pre[i] * x % P
    acc = pow(pre[n], P - 2, P)
    out = [0] * n
    for i in range(n - 1, -1, -1):
        out[i] = acc * pre[i] % P
        acc = acc * xs[i] % P
    return out

def batch_add(Ps, Qs):
    """affine add, distinct x assumed; returns list (None if exceptional)."""
    dens, idx = [], []
    out = [None] * len(Ps)
    for t, (a, b) in enumerate(zip(Ps, Qs)):
        if a is None: out[t] = b; continue
        if b is None: out[t] = a; continue
        if a[0] == b[0]:
            out[t] = ('EXC', a, b); continue
        dens.append((b[0] - a[0]) % P); idx.append(t)
    iv = batch_inv(dens) if dens else []
    for t, d in zip(idx, iv):
        ax, ay = Ps[t]; bx, by = Qs[t]
        l = (by - ay) * d % P
        cx = (l * l - ax - bx) % P
        out[t] = (cx, (l * (ax - cx) - ay) % P)
    return out
