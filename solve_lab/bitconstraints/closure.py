#!/usr/bin/env python3
"""Pass 4: forward closure from the 256 selector bits.

Given a concrete 0/1 assignment to the 256 selector bits, propagate through the
atom set (every atom is taken as a residual that must vanish -- the standing
model from METHOD_SUMMARY.md, re-validated here by the fact that the propagation
is contradiction-free and determines ~all wires).  Then evaluate every one of
the 39033 equations *exactly* and report which ones are violated.

An equation that is violated for some bit vectors and satisfied for others is a
genuine constraint on the bits.  An equation that is satisfied for *every* bit
vector carries no information about the bits.
"""
import pickle, json, os, sys, random, time
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
NVARS = 38748


def load():
    with open(os.path.join(HERE, 'cache.pkl'), 'rb') as f:
        D = pickle.load(f)
    chain = json.load(open(os.path.join(HERE, '..', 'anneal', 'chain.json')))['chain_bit_vars']
    return D, chain


def atom_vars(poly):
    s = set()
    for m, c in poly:
        s.update(m)
    return s


def reduce_poly(poly, val):
    """substitute known values; return dict monomial(unknown part)->coef"""
    out = defaultdict(int)
    for m, c in poly:
        rest = []
        coef = c
        for v in m:
            x = val[v]
            if x is None:
                rest.append(v)
            else:
                coef *= x
                if coef == 0:
                    break
        if coef == 0:
            continue
        out[tuple(rest)] += coef
    return {m: c for m, c in out.items() if c != 0}


def propagate(atoms, val, verbose=True):
    """Returns (n_determined, contradictions, undetermined_atoms)."""
    watch = defaultdict(list)
    for i, a in enumerate(atoms):
        for v in atom_vars(a):
            watch[v].append(i)
    queue = list(range(len(atoms)))
    inq = [True] * len(atoms)
    contradictions = []
    rounds = 0
    while queue:
        rounds += 1
        i = queue.pop()
        inq[i] = False
        red = reduce_poly(atoms[i], val)
        uv = set()
        for m in red:
            uv.update(m)
        if not uv:
            if red:   # nonzero constant
                contradictions.append(i)
            continue
        if len(uv) != 1:
            continue
        u = next(iter(uv))
        c = [0, 0, 0, 0]
        deg_ok = True
        for m, k in red.items():
            d = len(m)
            if d > 3:
                deg_ok = False
                break
            c[d] += k
        if not deg_ok:
            continue
        newval = None
        if c[2] == 0 and c[3] == 0:
            if c[1] == 0:
                if c[0] != 0:
                    contradictions.append(i)
                continue
            if c[0] % c[1] != 0:
                contradictions.append(i)
                continue
            newval = -c[0] // c[1]
        elif c[3] == 0:
            # quadratic c2 u^2 + c1 u + c0 = 0
            disc = c[1] * c[1] - 4 * c[2] * c[0]
            if disc < 0:
                contradictions.append(i)
                continue
            r = int(disc ** 0.5)
            while r * r > disc:
                r -= 1
            while (r + 1) * (r + 1) <= disc:
                r += 1
            if r * r != disc:
                contradictions.append(i)
                continue
            roots = set()
            for s in (r, -r):
                num = -c[1] + s
                if num % (2 * c[2]) == 0:
                    roots.add(num // (2 * c[2]))
            if len(roots) == 1:
                newval = roots.pop()
            else:
                continue
        else:
            continue
        if newval is None:
            continue
        val[u] = newval
        for j in watch[u]:
            if not inq[j]:
                inq[j] = True
                queue.append(j)
    det = sum(1 for v in val if v is not None)
    if verbose:
        print(f"  propagation: {rounds} atom visits, {det}/{NVARS} vars determined, "
              f"{len(contradictions)} contradicting atoms")
    return det, contradictions


def eval_eq(poly, val):
    """Exact value, or None if any needed var is unknown."""
    tot = 0
    for m, c in poly:
        t = c
        for v in m:
            x = val[v]
            if x is None:
                return None
            t *= x
            if t == 0:
                break
        tot += t
    return tot


def run(bitvec, D, chain, label):
    atoms = D['atoms']
    val = [None] * NVARS
    for b, v in zip(chain, bitvec):
        val[b] = v
    t0 = time.time()
    det, contra = propagate(atoms, val)
    viol, unknown = [], []
    for i, p in enumerate(D['eq_poly']):
        r = eval_eq(p, val)
        if r is None:
            unknown.append(i)
        elif r != 0:
            viol.append(i)
    print(f"[{label}] det={det} contra_atoms={len(contra)} "
          f"violated_eqs={len(viol)} undecided_eqs={len(unknown)} "
          f"({time.time()-t0:.1f}s)")
    return val, det, contra, viol, unknown


def main():
    D, chain = load()
    seeds = [('zeros', [0] * 256), ('ones', [1] * 256)]
    rng = random.Random(12345)
    for k in range(3):
        seeds.append((f'rand{k}', [rng.randrange(2) for _ in range(256)]))
    seeds.append(('e0', [1] + [0] * 255))
    results = {}
    for label, bv in seeds:
        val, det, contra, viol, unknown = run(bv, D, chain, label)
        results[label] = {'bits': bv, 'det': det,
                          'contra_atoms': contra[:200], 'n_contra': len(contra),
                          'violated': viol[:400], 'n_violated': len(viol),
                          'undecided': unknown[:400], 'n_undecided': len(unknown)}
    json.dump(results, open(os.path.join(HERE, 'closure_results.json'), 'w'))
    print("wrote closure_results.json")


if __name__ == '__main__':
    main()
