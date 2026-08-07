#!/usr/bin/env python3
"""EXACT verification of the confined-knob census.

jpins2.py used first-order derivatives.  A constraint with zero derivative can still
be disturbed at second order (atoms are degree 2), so the derivative census is only a
CANDIDATE list.  Here each candidate is verified by an actual finite perturbation:
release the pin, set the variable to several genuinely different values, re-propagate
exactly, and evaluate every constraint.  Confined only if NO constraint outside the
residual moves, at every probe.
"""
import sys, os, pickle, random, time
from collections import defaultdict, deque
import jengine as E, jman as J, jmodp as MP, jsolve2 as S
import jdist as DI

P = MP.P
HERE = os.path.dirname(os.path.abspath(__file__))
definer, order = J.definer, J.order
polys = E.polys
pos = {v: k for k, v in enumerate(order)}
FREE = set(range(E.NV)) - set(definer)
uses = defaultdict(list)
for v, i in definer.items():
    for w in E.varsof[i]:
        if w != v:
            uses[w].append(v)
CONS = sorted(set(range(len(polys))) - set(definer.values()))
EVP = {}
for v, i in definer.items():
    p = polys[i]
    EVP[v] = (pow(p[(v,)] % P, P - 2, P),
              tuple((k, cc % P) for k, cc in p.items() if k != (v,)))
BAD = [20407, 20409, 31575]
BADSET = set(BAD)
base, _ = S.branch(1, 1)
R0 = {i: MP.atom_modp(i, base) for i in CONS}


def cone_of(v):
    c = set(); q = deque([v])
    while q:
        x = q.popleft()
        for w in uses[x]:
            if w != v and w not in c and w not in FREE:
                c.add(w); q.append(w)
    return sorted(c, key=lambda z: pos[z])


def exact_disturbed(v, deltas):
    """actual finite perturbation of a RELEASED v; returns constraints that move"""
    cone = cone_of(v)
    touched = set()
    for w in [v] + cone:
        touched.update(i for i in CONS if w in E.varsof[i])
    moved = set()
    for d in deltas:
        val = list(base)
        val[v] = (val[v] + d) % P
        for w in cone:
            e = EVP.get(w)
            if e is None:
                continue
            ic, rest = e
            s = 0
            for k, cc in rest:
                t = cc
                for j in k:
                    t = t * val[j] % P
                s += t
            val[w] = (-s) % P * ic % P
        for i in touched:
            if MP.atom_modp(i, val) != R0[i]:
                moved.add(i)
    return moved


if __name__ == '__main__':
    cand = pickle.load(open(os.path.join(HERE, 'jpins2.pkl'), 'rb'))
    cand.sort()
    random.seed(17)
    deltas = [1, 2, 7, random.randrange(P), random.randrange(P)]
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 40
    print(f"exact-verifying the {min(n,len(cand))} cheapest of {len(cand)} "
          f"first-order-confined candidates")
    survivors = []
    t0 = time.time()
    for cost, v, pin, firstorder in cand[:n]:
        mv = exact_disturbed(v, deltas)
        ok = mv <= BADSET and mv
        tag = "CONFINED" if ok else f"leaks {len(mv - BADSET)}"
        print(f"  cost {cost:3d} x_{v:<7} pin a{pin:<7} first-order {list(firstorder)} "
              f"-> exact {sorted(mv & BADSET)}  {tag}", flush=True)
        if ok:
            survivors.append((cost, v, pin, tuple(sorted(mv))))
    print(f"\nEXACT confined knobs among the {n} cheapest: {len(survivors)}  "
          f"({time.time()-t0:.0f}s)")
    for s in survivors:
        print("   ", s)
    pickle.dump(survivors, open(os.path.join(HERE, 'jpins3.pkl'), 'wb'))
