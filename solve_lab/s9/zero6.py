"""Only the two setter pins remain on the all-zero branch:
    688  : 8863713*(x_18956 - C2const) = x_14257 = p*handle   ->  x_18956 = C2const (mod p)
    1618 : x_24468 - C1const = x_32989 = p*handle             ->  x_24468 = C1const (mod p)
Both carriers are movable.  Measure the sensitivities and solve, then set the handles in Z.
"""
import pickle, sys
import harness as H
exec(open('repair.py').read().split('if __name__')[0])
P = 2**256 - 2**32 - 977
d = pickle.load(open('atoms.pkl', 'rb')); src = d['atom_src']
a2e = {a: set(e) for a, e in pickle.load(open('atom2eq.pkl', 'rb')).items()}
C1c = 91416258160755509149180373473728639746431157665678710450404458852172057265575180278101002
C2c = 125787314747601108116039725163361763116550465675981151838811516827327919228823597744635626
NV = 38748
freeinp = [x for x in range(NV) if x not in definer]
CODES, _ = H.load_equations()


def allnz(v):
    return sorted(a for a, Pp in enumerate(polys) if evalpoly(Pp, v) != 0)


if __name__ == '__main__':
    v0 = H.load_assignment('zero5_out.json')
    print('start atoms:', allnz(v0), 'failing =', len(H.evaluate(CODES, v0)))
    print(f'x_18956 %p = {v0[18956] % P}\n   target C2c %p = {C2c % P}')
    print(f'x_24468 %p = {v0[24468] % P}\n   target C1c %p = {C1c % P}')
    # which free inputs move each carrier, and by how much per unit?
    mv = {18956: [], 24468: []}
    for f in freeinp:
        w = list(v0); ripple(w, {f: v0[f] + 1})
        for t in (18956, 24468):
            if w[t] != v0[t]:
                mv[t].append((f, (w[t] - v0[t]) % P))
    for t in (18956, 24468):
        print(f'\nfree inputs moving x_{t}: {len(mv[t])}')
        for f, dd in mv[t][:8]: print(f'   x_{f}: d={dd}')
    # try single-knob solves
    print()
    for t, target in ((18956, C2c), (24468, C1c)):
        need = (target - v0[t]) % P
        for f, dd in mv[t]:
            if dd == 0: continue
            k = need * pow(dd, P - 2, P) % P
            w = list(v0); ripple(w, {f: v0[f] + k})
            if (w[t] - target) % P == 0:
                nz = allnz(w)
                print(f'x_{t}: knob x_{f} with k -> residue MATCHED; atoms now {nz}')
                break
        else:
            print(f'x_{t}: no single knob reaches the target residue')
