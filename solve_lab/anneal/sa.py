#!/usr/bin/env python3
"""sa.py -- a plain simulated annealer over a QUBO, used as a classical stand-in for the QPU."""
import random, math


class SA:
    def __init__(self, Q, nvars):
        self.n = nvars
        self.h = [0.0] * nvars
        self.adj = [[] for _ in range(nvars)]
        self.off = 0.0
        for m, c in Q.items():
            if not m: self.off += c
            elif len(m) == 1: self.h[m[0]] += c
            else:
                i, j = m
                self.adj[i].append((j, c)); self.adj[j].append((i, c))

    def energy(self, x):
        e = self.off + sum(self.h[i] for i in range(self.n) if x[i])
        for i in range(self.n):
            if x[i]:
                for j, c in self.adj[i]:
                    if j > i and x[j]: e += c
        return e

    def run(self, sweeps=2000, beta0=0.05, beta1=30.0, seed=0, x0=None, clamp=()):
        rnd = random.Random(seed)
        x = list(x0) if x0 else [rnd.randrange(2) for _ in range(self.n)]
        clamped = set(clamp)
        free = [i for i in range(self.n) if i not in clamped]
        d = [self.h[i] + sum(c for j, c in self.adj[i] if x[j]) for i in range(self.n)]
        e = self.energy(x)
        best, bx = e, list(x)
        for s in range(sweeps):
            beta = beta0 * (beta1 / beta0) ** (s / max(1, sweeps - 1))
            rnd.shuffle(free)
            for i in free:
                de = d[i] if not x[i] else -d[i]
                if de <= 0 or rnd.random() < math.exp(-beta * de):
                    sgn = 1 if not x[i] else -1
                    x[i] ^= 1; e += de
                    for j, c in self.adj[i]: d[j] += sgn * c
                    if e < best: best, bx = e, list(x)
            if best == 0: break
        return best, bx
