#!/usr/bin/env python3
"""Agent Z: (i) global separability check -- does ANY monomial anywhere in the
instance contain two DISTINCT selector variables?  (ii) checker-anchored test of
the selector-only equations at every weight from 0 to 256."""
import os, sys, json, pickle, collections, random

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sel = set(json.load(open(os.path.join(HERE, 'zsel.json')))['selectors'])
eqp, boolvars = pickle.load(open(os.path.join(HERE, 'zbool.pkl'), 'rb'))

# ---- (i) separability
worst = 0
bad = []
seldeg = collections.Counter()
for i, p in enumerate(eqp):
    for m in p:
        s = set(m) & sel
        seldeg[len(s)] += 1
        if len(s) >= 2:
            bad.append((i, m))
        worst = max(worst, len(s))
print("monomials by number of DISTINCT selectors they contain:", sorted(seldeg.items()))
print("monomials containing >=2 distinct selectors:", len(bad))
print("=> every equation is AFFINE in the selector vector:", worst <= 1)

# selector total degree inside a monomial (s^2 shows up)
tot = collections.Counter()
for i, p in enumerate(eqp):
    for m in p:
        tot[sum(1 for v in m if v in sel)] += 1
print("selector total-degree per monomial (with multiplicity):", sorted(tot.items()))

# ---- equations that vanish identically once every boolean var is boolean
def red(p):
    q = {}
    for m, c in p.items():
        mm = tuple(sorted(v for v in set(m) if True) for _ in [0])[0]
        pass
    return q

ident = []
for i, p in enumerate(eqp):
    q = {}
    for m, c in p.items():
        mm = tuple(sorted(set(v for v in m) if all(v in boolvars for v in m) else m))
        # x^2 -> x only for boolean vars
        cnt = collections.Counter(m)
        mm = []
        for v, k in cnt.items():
            mm.extend([v] if v in boolvars else [v] * k)
        mm = tuple(sorted(mm))
        q[mm] = q.get(mm, 0) + c
    q = {m: c for m, c in q.items() if c}
    if not q:
        ident.append(i)
print("equations that vanish IDENTICALLY on the boolean locus:", len(ident))

# ---- (ii) checker-anchored: evaluate those equations under many sigma
sys.path.insert(0, os.path.join(HERE, '..'))
import importlib.util
spec = importlib.util.spec_from_file_location("chk", os.path.join(HERE, '..', 'checker.py'))
chk = importlib.util.module_from_spec(spec)
spec.loader.exec_module(chk)
codes, varsets = chk.load_equations()
print("checker loaded", len(codes), "equations")

base = chk.load_assignment(os.path.join(HERE, '..', 'best', 'new_instance_partial_39026.json'))
SELL = sorted(sel)
rng = random.Random(7)
print()
print("%5s  %s" % ("w", "satisfied among the %d identically-vanishing eqs (checker eval)" % len(ident)))
for w in [0, 1, 2, 3, 8, 17, 32, 64, 128, 192, 250, 255, 256]:
    on = set(rng.sample(SELL, w))
    v = list(base)
    for s in SELL:
        v[s] = 1 if s in on else 0
    ns = {'v': v, '__builtins__': {}}
    ok = sum(1 for i in ident if eval(codes[i], ns) == 0)
    # also: full-instance score under this brutal overwrite, for context
    print("%5d  %6d / %d" % (w, ok, len(ident)))

# the 13 pure-selector equations, at every weight 0..256 (random sigma each)
pure = []
for i, p in enumerate(eqp):
    vs = set()
    for m in p:
        vs |= set(m)
    if vs and vs <= sel:
        pure.append(i)
print()
print("pure-selector equations:", len(pure), pure)
allok = True
for w in range(0, 257):
    on = set(rng.sample(SELL, w))
    v = [0] * chk.NVARS
    for s in on:
        v[s] = 1
    ns = {'v': v, '__builtins__': {}}
    for i in pure:
        if eval(codes[i], ns) != 0:
            allok = False
            print("  VIOLATED at w=%d eq %d" % (w, i))
print("all 13 pure-selector equations satisfied at every weight 0..256 (checker eval):", allok)
