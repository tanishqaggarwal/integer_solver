"""Full all-zero-branch pipeline, parameterised by which gate-bit satisfies the forced OR gate.

The mirror chain from the setter pins 688/1618 terminates on x_24221 and x_25477.  Those are
pinned only when gate-bit x_47 is ON.  So choosing a DIFFERENT gate-bit for the OR gate should
let the chain terminate freely.  Try every gate-bit.
"""
import pickle, sys, time
import harness as H
exec(open('repair.py').read().split('if __name__')[0])
P = 2**256 - 2**32 - 977
d = pickle.load(open('atoms.pkl', 'rb')); src = d['atom_src']
a2e = {a: set(e) for a, e in pickle.load(open('atom2eq.pkl', 'rb')).items()}
boolv = set(pickle.load(open('boolvars.pkl', 'rb')))
pins = pickle.load(open('pins/pins.pkl', 'rb'))
bypin = {}
for pn in pins: bypin.setdefault(pn['G'], []).append(pn)
GATEBITS = sorted(bypin)
NV = 38748
freeinp = [x for x in range(NV) if x not in definer]
nbfree = [f for f in freeinp if f not in boolv]
CODES, _ = H.load_equations()
C1c = 91416258160755509149180373473728639746431157665678710450404458852172057265575180278101002
C2c = 125787314747601108116039725163361763116550465675981151838811516827327919228823597744635626


def allnz(v):
    return sorted(a for a, Pp in enumerate(polys) if evalpoly(Pp, v) != 0)


def one_driver(v, t):
    """the non-boolean free input driving x_t exactly 1:1 mod p, if unique-ish"""
    for f in nbfree:
        w = list(v); ripple(w, {f: v[f] + 1})
        if w[t] != v[t] and (w[t] - v[t]) % P == 1:
            return f
    return None


def atom_pattern(a):
    Pp = polys[a]
    lin = {m[0]: c for m, c in Pp.items() if len(m) == 1}
    if len(lin) != 3: return None
    items = sorted(lin.items(), key=lambda kv: -abs(kv[1]))
    (v1, c1), (v2, c2), (v3, c3) = items
    if c1 != -c2 or abs(c3) != 1: return None
    return (abs(c1), v1, v2, v3)


def set_handle(v, Hv, num):
    hd = definer.get(Hv)
    if hd is None or num % P: return
    hp = polys[hd]
    for m in hp:
        if len(m) == 2 and Hv not in m:
            w1, w2 = m
            base = w1 if v[w2] == P else (w2 if v[w1] == P else None)
            if base is not None:
                ripple(v, {base: num // P}); return


def pipeline(G, maxchase=25, verbose=False):
    v = H.load_assignment('zero_out.json')
    ripple(v, {22162: 0, 30213: 0})
    seeds = {G: 1}
    for pn in bypin[G]: seeds[pn['B']] = pn['HUGE']
    ripple(v, seeds)
    # close setter pins 688 / 1618
    for carrier, target, coef, Hv in ((18956, C2c, 8863713, 14257), (24468, C1c, 1, 32989)):
        drv = one_driver(v, carrier)
        if drv is None: continue
        ripple(v, {drv: v[drv] + (target - v[carrier]) % P})
        set_handle(v, Hv, coef * (v[carrier] - target))
    # chase the mirror chain
    seen = set()
    for _ in range(maxchase):
        nz = allnz(v)
        if not nz: break
        prim = [a for a in nz if atom_pattern(a)]
        moved = False
        for a in prim:
            coef, X, Y, Hv = atom_pattern(a)
            for t, other in ((Y, X), (X, Y)):
                drv = one_driver(v, t)
                if drv is None or (a, drv) in seen: continue
                seen.add((a, drv))
                ripple(v, {drv: v[drv] + (v[other] - v[t]) % P})
                set_handle(v, Hv, coef * (v[X] - v[Y]))
                moved = True; break
            if moved: break
        if not moved: break
    return v


if __name__ == '__main__':
    sel = [int(x) for x in sys.argv[1:]] or GATEBITS
    best = None; t0 = time.time()
    for i, G in enumerate(sel):
        try:
            v = pipeline(G)
        except Exception as e:
            print(f'x_{G}: ERROR {e}'); continue
        nz = allnz(v); f = H.evaluate(CODES, v)
        print(f'x_{G}: {len(f)} failing, atoms {nz}', flush=True)
        if best is None or len(f) < best[0]:
            best = (len(f), G, nz)
            H.save_assignment(v, 'branch_best.json')
            if not f:
                print('*** ZERO FAILING EQUATIONS ***'); break
    print('\nBEST:', best)
