"""Close the last two setter pins on the all-zero branch.
x_16742 moves x_18956 by exactly +1 mod p; x_14681 moves x_24468 by exactly +1 mod p.
Both are non-boolean free inputs.  Set them to hit the pinned residues, then set the
quotient handles so the pins hold exactly in Z."""
import pickle, sys
import harness as H
exec(open('repair.py').read().split('if __name__')[0])
P = 2**256 - 2**32 - 977
d = pickle.load(open('atoms.pkl', 'rb')); src = d['atom_src']
a2e = {a: set(e) for a, e in pickle.load(open('atom2eq.pkl', 'rb')).items()}
C1c = 91416258160755509149180373473728639746431157665678710450404458852172057265575180278101002
C2c = 125787314747601108116039725163361763116550465675981151838811516827327919228823597744635626
CODES, _ = H.load_equations()


def allnz(v):
    return sorted(a for a, Pp in enumerate(polys) if evalpoly(Pp, v) != 0)


def stage(v, tag):
    nz = allnz(v); f = H.evaluate(CODES, v)
    print(f'[{tag}] atoms={nz}  EQ {len(CODES)-len(f)}/{len(CODES)} ({len(f)} failing)')
    return nz, f


if __name__ == '__main__':
    v = H.load_assignment('zero5_out.json')
    stage(v, '0 start')
    d1 = (C2c - v[18956]) % P
    d2 = (C1c - v[24468]) % P
    print(f'shift x_16742 by {d1}\nshift x_14681 by {d2}')
    ripple(v, {16742: v[16742] + d1, 14681: v[14681] + d2})
    print(f'  x_18956 = C2c mod p ? {(v[18956]-C2c) % P == 0}')
    print(f'  x_24468 = C1c mod p ? {(v[24468]-C1c) % P == 0}')
    stage(v, '1 residues matched')
    # now the handles: 688 needs x_14257 = 8863713*(x_18956-C2c); 1618 needs x_32989 = x_24468-C1c
    n688 = 8863713 * (v[18956] - C2c)
    n1618 = v[24468] - C1c
    print(f'  688 handle numerator divisible by p ? {n688 % P == 0}')
    print(f'  1618 handle numerator divisible by p ? {n1618 % P == 0}')
    seeds = {}
    if n688 % P == 0: seeds[7497] = n688 // P
    if n1618 % P == 0: seeds[11436] = n1618 // P
    ripple(v, seeds)
    stage(v, '2 handles set')
    ok, hist = repair_loop(v, rounds=20, verbose=False)
    nz, f = stage(v, '3 after repair')
    for a in nz:
        print(f'    atom {a} ({len(a2e.get(a,[]))} eqs) gate={atom_out.get(a)}: {src[a][:95]}')
    H.save_assignment(v, 'zero7_out.json')
