#!/usr/bin/env python3
"""
SOLVE the concentrated core instead of annealing it.

concentrate2.py showed the hardness concentrates into

    A * c^2  ==  B^2   (mod p)

with A, B, c AFFINE (measured, 4/4 additive) in four knobs
x14853, x16742, x22162, x22649.  Affine inputs + a quadratic condition means the
core is a LOW-DEGREE POLYNOMIAL over GF(p): fix three knobs, and the fourth
satisfies a univariate cubic.  Univariate roots mod p are found in polynomial
time (Cantor-Zassenhaus), so this block does not need an annealer at all --
which is the real payoff of concentrating the hardness here.

Steps: fit the affine models -> build the univariate polynomial -> find its roots
mod p -> realise a root -> EXACT integer check.

Usage: python3 solve_core.py
"""
import os, sys, time, json, random

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 's9', 'eff'))
sys.path.insert(0, os.path.join(HERE, '..', 's10'))
sys.path.insert(0, HERE)
import lib as L
import ad

LAB = os.path.join(HERE, '..')
P = 2**256 - 2**32 - 977
CORE = {'A': 33469, 'B': 27713, 'c': 1326, 'u': 29322, 'w': 3558}
KNOBS = [14853, 16742, 22162, 22649]


# ------------------------------------------------ polynomials mod p ---------
def pmul(f, g):
    r = [0] * (len(f) + len(g) - 1)
    for i, a in enumerate(f):
        if a:
            for j, b in enumerate(g):
                r[i + j] = (r[i + j] + a * b) % P
    return trim(r)


def padd(f, g):
    n = max(len(f), len(g))
    return trim([( (f[i] if i < len(f) else 0) + (g[i] if i < len(g) else 0)) % P
                 for i in range(n)])


def psub(f, g):
    n = max(len(f), len(g))
    return trim([((f[i] if i < len(f) else 0) - (g[i] if i < len(g) else 0)) % P
                 for i in range(n)])


def trim(f):
    while len(f) > 1 and f[-1] % P == 0:
        f.pop()
    return f


def pmod(f, g):
    f = f[:]
    dg = len(g) - 1
    inv = pow(g[-1], -1, P)
    while len(f) - 1 >= dg and any(f):
        d = len(f) - 1 - dg
        co = f[-1] * inv % P
        for i in range(len(g)):
            f[d + i] = (f[d + i] - co * g[i]) % P
        trim(f)
        if len(f) - 1 < dg:
            break
    return trim(f)


def pgcd(f, g):
    while len(g) > 1 or g[0] % P:
        f, g = g, pmod(f, g)
        if len(g) == 1 and g[0] % P == 0:
            break
    return f


def ppowmod(base, e, mod):
    result = [1]
    b = pmod(base[:], mod)
    while e:
        if e & 1:
            result = pmod(pmul(result, b), mod)
        b = pmod(pmul(b, b), mod)
        e >>= 1
    return result


def roots_mod_p(f):
    """All roots of f in GF(p) (f of small degree)."""
    f = trim([c % P for c in f])
    if len(f) == 1:
        return []
    # gcd(x^p - x, f) = product of distinct linear factors
    xp = ppowmod([0, 1], P, f)
    g = pgcd(f[:], psub(xp, [0, 1]))
    if len(g) == 1:
        return []
    out = []

    def split(h):
        h = trim(h)
        if len(h) == 1:
            return
        if len(h) == 2:                       # a*x + b
            out.append((-h[0] * pow(h[1], -1, P)) % P)
            return
        for _ in range(60):
            a = random.randrange(P)
            t = ppowmod([a, 1], (P - 1) // 2, h)
            d = pgcd(h[:], psub(t, [1]))
            if 1 < len(d) < len(h):
                split(d); split(pmod(h, d) if False else pdiv(h, d))
                return
        return

    def pdiv(f, g):
        f = f[:]; q = [0] * (len(f) - len(g) + 1)
        inv = pow(g[-1], -1, P)
        while len(f) - 1 >= len(g) - 1 and any(f):
            d = len(f) - 1 - (len(g) - 1)
            co = f[-1] * inv % P
            q[d] = co
            for i in range(len(g)):
                f[d + i] = (f[d + i] - co * g[i]) % P
            trim(f)
            if len(f) - 1 < len(g) - 1:
                break
        return trim(q)

    split(g)
    return sorted(set(out))


# ------------------------------------------------------------- main ---------
def main():
    t0 = time.time()
    random.seed(7)
    v0 = L.load(os.path.join(LAB, 'best', 'new_instance_partial_39026.json'))

    def core_at(delta):
        vv = list(v0)
        for x, d in delta.items():
            vv[x] += d
        ad.fwd(vv, rounds=6)
        return {k: vv[i] % P for k, i in CORE.items()}

    base = core_at({})
    print("affine models for A, B, c in the four knobs:")
    coef = {}
    for x in KNOBS:
        e = core_at({x: 1})
        coef[x] = {k: (e[k] - base[k]) % P for k in ('A', 'B', 'c')}
        nz = {k: ('0' if coef[x][k] == 0 else 'nonzero') for k in ('A', 'B', 'c')}
        print(f"  d/d x{x}: A {nz['A']}, B {nz['B']}, c {nz['c']}")

    # verify the affine model before trusting it
    d = {x: random.randrange(1, 10**9) for x in KNOBS}
    act = core_at(d)
    okmodel = True
    for k in ('A', 'B', 'c'):
        pred = (base[k] + sum(coef[x][k] * d[x] for x in KNOBS)) % P
        if pred != act[k]:
            okmodel = False
    print(f"  affine model verified at a random point: {okmodel}")
    if not okmodel:
        print("  model invalid -> stopping (would be meaningless)")
        return

    # ---- univariate: vary one knob, others fixed ---------------------------
    print("\nsolving  A*c^2 - B^2 == 0 (mod p)  as a univariate polynomial")
    for x in KNOBS:
        aA = [base['A'], coef[x]['A']]
        aB = [base['B'], coef[x]['B']]
        aC = [base['c'], coef[x]['c']]
        f = psub(pmul(aA, pmul(aC, aC)), pmul(aB, aB))
        deg = len(f) - 1
        r = roots_mod_p(f[:]) if deg >= 1 else []
        print(f"  knob x{x}: polynomial degree {deg}, roots mod p: {len(r)}")
        for t in r[:2]:
            # verify the root really zeroes the core
            chk = core_at({x: t})
            disc = (chk['A'] * chk['c'] * chk['c'] - chk['B'] * chk['B']) % P
            print(f"     root {str(t)[:30]}... -> A c^2 - B^2 = "
                  f"{'0  CORE SOLVED' if disc == 0 else 'nonzero (model/realisation gap)'}")
            if disc == 0:
                vv = list(v0); vv[x] += t; ad.fwd(vv, rounds=6)
                av = L.all_atom_values(vv)
                sc = L.NEQ - len(L.failing_eqs(av))
                print(f"     realised score: {sc}")
                out = os.path.join(HERE, f'core_root_x{x}.json')
                json.dump({f'x_{i}': int(vv[i]) for i in range(len(vv)) if vv[i]},
                          open(out, 'w'))
                print(f"     written -> {out}")
    print(f"\n{time.time()-t0:.0f}s")


if __name__ == '__main__':
    main()
