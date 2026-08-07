#!/usr/bin/env python3
"""Full confined-knob census.

Breaking the definer atom of v costs |eqs(definer(v))| equations and frees v.
The knob is USABLE only if CONFINED: perturbing v (with the whole forward cone
re-derived) must disturb no constraint outside the residual.

Releasing a pin == seeding der[v]=1 and refusing to let v's own definer overwrite it,
then propagating derivatives forward.  So one RESTRICTED forward-AD pass per candidate
(only v's forward cone) answers it; no DAG rebuild needed.

Scans every defined variable that moves the residual at all (from jrev.pkl).
"""
import sys, os, pickle, time
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

# constraint atoms, and which constraints each variable can appear in
CONS = sorted(set(range(len(polys))) - set(definer.values()))
CONSSET = set(CONS)
var2cons = defaultdict(list)
for i in CONS:
    for v in E.varsof[i]:
        var2cons[v].append(i)

EVP = {}
for v, i in definer.items():
    p = polys[i]
    EVP[v] = (pow(p[(v,)] % P, P - 2, P),
              tuple((k, cc % P) for k, cc in p.items() if k != (v,)))

base, _ = S.branch(1, 1)
BAD = [20407, 20409, 31575]
BADSET = set(BAD)


def datom(i, val, der):
    s = 0
    for k, c in polys[i].items():
        t = c % P
        dt = 0
        for j in k:
            dt = (dt * val[j] + t * der.get(j, 0)) % P
            t = t * val[j] % P
        s += dt
    return s % P


def disturbed(v, val):
    """constraints with nonzero gradient when v is released and perturbed"""
    cone = set()
    q = deque([v])
    while q:
        x = q.popleft()
        for w in uses[x]:
            if w != v and w not in cone and w not in FREE:
                cone.add(w); q.append(w)
    der = {v: 1}
    for w in sorted(cone, key=lambda z: pos[z]):
        e = EVP.get(w)
        if e is None:
            continue
        ic, rest = e
        ds = 0
        for k, cc in rest:
            t = cc
            dt = 0
            for j in k:
                dt = (dt * val[j] + t * der.get(j, 0)) % P
                t = t * val[j] % P
            ds += dt
        d = (-ds) % P * ic % P
        if d:
            der[w] = d
    touched = set()
    for w in list(der):
        touched.update(var2cons.get(w, ()))
    out = []
    for i in touched:
        if datom(i, val, der):
            out.append(i)
    return out


if __name__ == '__main__':
    R = pickle.load(open(os.path.join(HERE, 'jrev.pkl'), 'rb'))
    grads = R['grads']
    movers = set()
    for c in grads:
        movers |= set(grads[c])
    cands = [(len(DI.A2E[definer[v]]), v, definer[v]) for v in movers if v in definer]
    cands.sort()
    print(f"defined movers (breakable pins) to scan: {len(cands)}")

    t0 = time.time()
    confined = []
    for n, (cost, v, pin) in enumerate(cands):
        d = disturbed(v, base)
        if set(d) <= BADSET and d:
            confined.append((cost, v, pin, tuple(sorted(d))))
            print(f"  CONFINED: cost {cost:3d} x_{v} pin a{pin} moves {sorted(d)}", flush=True)
        if n % 200 == 0:
            print(f"   ... {n}/{len(cands)} scanned, {len(confined)} confined "
                  f"({time.time()-t0:.0f}s)", flush=True)
    print(f"\nSCAN COMPLETE: {len(cands)} pins scanned, {len(confined)} CONFINED")
    for c in confined:
        print("   ", c)
    pickle.dump(confined, open(os.path.join(HERE, 'jpins2.pkl'), 'wb'))
