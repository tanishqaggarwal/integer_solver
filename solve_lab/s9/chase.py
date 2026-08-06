"""Automated mirror-chase on the all-zero branch.

Every remaining obstruction has the shape   k*(F - C) - handle = 0   where F is a free input we
just moved and C is a computed mirror.  Close it by finding the non-boolean free input that drives
C exactly 1:1 mod p, moving it by the required delta, and setting the handle in Z.  That in turn
lights the mirror one level up.  Chase it and see whether the chain terminates or cycles.
"""
import pickle, sys, time
import harness as H
exec(open('repair.py').read().split('if __name__')[0])
P = 2**256 - 2**32 - 977
d = pickle.load(open('atoms.pkl', 'rb')); src = d['atom_src']
a2e = {a: set(e) for a, e in pickle.load(open('atom2eq.pkl', 'rb')).items()}
boolv = set(pickle.load(open('boolvars.pkl', 'rb')))
NV = 38748
freeinp = [x for x in range(NV) if x not in definer]
nbfree = [f for f in freeinp if f not in boolv]
CODES, _ = H.load_equations()


def allnz(v):
    return sorted(a for a, Pp in enumerate(polys) if evalpoly(Pp, v) != 0)


def drivers(v, targets):
    """non-boolean free inputs moving each target, with per-unit delta mod p"""
    out = {t: [] for t in targets}
    for f in nbfree:
        w = list(v); ripple(w, {f: v[f] + 1})
        for t in targets:
            if w[t] != v[t]:
                dd = (w[t] - v[t]) % P
                if dd: out[t].append((f, dd))
    return out


def atom_pattern(a):
    """for a check atom, return (coef, freeVar, mirrorVar, handleVar) if it matches k*(X-Y)-H"""
    Pp = polys[a]
    lin = {m[0]: c for m, c in Pp.items() if len(m) == 1}
    if len(lin) != 3: return None
    items = sorted(lin.items(), key=lambda kv: -abs(kv[1]))
    (v1, c1), (v2, c2), (v3, c3) = items
    if c1 != -c2: return None
    if abs(c3) != 1: return None
    return (abs(c1), v1, v2, v3)


if __name__ == '__main__':
    v = H.load_assignment(sys.argv[1] if len(sys.argv) > 1 else 'zero9_out.json')
    seen_pairs = set(); chain = []
    for step in range(40):
        nz = allnz(v)
        f = H.evaluate(CODES, v)
        print(f'--- step {step}: atoms={nz}  failing={len(f)}', flush=True)
        if not nz:
            print('*** ALL ATOMS ZERO ***'); break
        prim = [a for a in nz if atom_pattern(a)]
        if not prim:
            print('   no mirror-shaped atom left; stopping'); break
        a = prim[0]
        pat = atom_pattern(a)
        coef, X, Y, Hv = pat
        # orient: the side that a non-boolean free input drives 1:1 is the mirror we move
        cand = drivers(v, [X, Y])
        pick = None
        for t in (Y, X):
            ones = [(fv, dd) for fv, dd in cand[t] if dd == 1]
            if ones: pick = (t, ones[0][0]); break
        if pick is None:
            print(f'   atom {a}: no 1:1 driver for {X} or {Y}; stopping'); break
        t, fv = pick
        other = X if t == Y else Y
        key = (a, fv)
        if key in seen_pairs:
            print(f'   CYCLE detected at atom {a} via x_{fv}; stopping'); break
        seen_pairs.add(key)
        delta = (v[other] - v[t]) % P
        chain.append((a, fv, t, other))
        print(f'   atom {a}: move x_{fv} to drive x_{t} -> x_{other} (delta {str(delta)[:24]}...)', flush=True)
        ripple(v, {fv: v[fv] + delta})
        # set the handle exactly
        num = coef * (v[X] - v[Y]) if abs(v[X] - v[Y]) else 0
        hd = definer.get(Hv)
        if hd is not None:
            hp = polys[hd]
            prods = [m for m in hp if len(m) == 2 and Hv not in m]
            if prods:
                m = prods[0]
                w1, w2 = m
                base = w1 if v[w2] == P else (w2 if v[w1] == P else None)
                if base is not None and num % P == 0:
                    ripple(v, {base: num // P})
    print('\nchain length:', len(chain))
    for a, fv, t, other in chain: print(f'   atom {a}: x_{fv} -> x_{t} (mirror of x_{other})')
    H.save_assignment(v, 'chase_out.json')
    print('final failing:', len(H.evaluate(CODES, v)))
