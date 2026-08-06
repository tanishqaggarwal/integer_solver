"""On the all-zero branch only the two setter pins 688/1618 and the forced OR gate remain.
The OR gate needs some boolean set, but 256 of the 1,156 booleans are GATE-BITS carrying two
load pins each -- flipping those cascades.  The other ~900 carry no pins at all.
Find an OR witness among the pin-free booleans."""
import pickle, sys
import harness as H
exec(open('repair.py').read().split('if __name__')[0])
P = 2**256 - 2**32 - 977
d = pickle.load(open('atoms.pkl', 'rb')); src = d['atom_src']
a2e = {a: set(e) for a, e in pickle.load(open('atom2eq.pkl', 'rb')).items()}
boolv = set(pickle.load(open('boolvars.pkl', 'rb')))
pins = pickle.load(open('pins/pins.pkl', 'rb'))
GATEBITS = set(pn['G'] for pn in pins)
NV = 38748
freeinp = [x for x in range(NV) if x not in definer]
bfree = [b for b in freeinp if b in boolv]
nopin = [b for b in bfree if b not in GATEBITS]


def allnz(v):
    return sorted(a for a, Pp in enumerate(polys) if evalpoly(Pp, v) != 0)


def score(v):
    codes, _ = H.load_equations()
    return len(H.evaluate(codes, v))


if __name__ == '__main__':
    print(f'boolean free inputs: {len(bfree)};  gate-bits (with pins): {len(GATEBITS & set(bfree))};'
          f'  PIN-FREE booleans: {len(nopin)}')
    v0 = H.load_assignment('zero_out.json')
    ripple(v0, {22162: 0, 30213: 0})
    print('base (all-zero branch, carriers zeroed):', allnz(v0), 'failing =', score(v0))
    results = []
    for b in nopin:
        w = list(v0); ripple(w, {b: 1 - v0[b]})
        if w[9274] != 1 or w[15298] != 0:
            continue
        nz = allnz(w)
        results.append((len(nz), b, nz))
    results.sort()
    print(f'pin-free booleans that satisfy the OR gate with the core off: {len(results)}')
    for n, b, nz in results[:10]:
        print(f'   x_{b}: {n} atoms {nz}')
    if results:
        n, b, nz = results[0]
        w = list(v0); ripple(w, {b: 1 - v0[b]})
        ok, hist = repair_loop(w, rounds=15, verbose=False)
        nz = allnz(w)
        f = score(w)
        print(f'\nbest pin-free OR witness x_{b}: atoms={nz}  failing={f}')
        for a in nz:
            print(f'    atom {a} ({len(a2e.get(a,[]))} eqs) gate={atom_out.get(a)}: {src[a][:95]}')
        H.save_assignment(w, 'zero4_out.json')
