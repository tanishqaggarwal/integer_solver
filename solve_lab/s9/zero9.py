"""Close the last two mirror pins on the all-zero branch.
    26731 : 6788513*(x_16742 - x_19083) = x_9254  = p*x_33787   -> move x_19083 via free x_38667
    9193  : 14572767*(x_14681 - x_24483) = x_18063 = ?*x_19890  -> move x_24483 via free x_29851
"""
import pickle, sys
import harness as H
exec(open('repair.py').read().split('if __name__')[0])
P = 2**256 - 2**32 - 977
d = pickle.load(open('atoms.pkl', 'rb')); src = d['atom_src']
a2e = {a: set(e) for a, e in pickle.load(open('atom2eq.pkl', 'rb')).items()}
CODES, _ = H.load_equations()


def allnz(v):
    return sorted(a for a, Pp in enumerate(polys) if evalpoly(Pp, v) != 0)


def stage(v, tag):
    nz = allnz(v); f = H.evaluate(CODES, v)
    print(f'[{tag}] atoms={nz}  EQ {len(CODES)-len(f)}/{len(CODES)} ({len(f)} failing)')
    return nz, f


if __name__ == '__main__':
    v = H.load_assignment('zero8_out.json')
    stage(v, '0 start')
    d1 = (v[16742] - v[19083]) % P          # x_19083 must rise by this to match x_16742
    d2 = (v[14681] - v[24483]) % P
    print(f'move x_19083 by {d1}\nmove x_24483 by {d2}')
    ripple(v, {38667: v[38667] + d1, 29851: v[29851] + d2})
    print(f'  x_16742 = x_19083 mod p ? {(v[16742]-v[19083]) % P == 0}')
    print(f'  x_14681 = x_24483 mod p ? {(v[14681]-v[24483]) % P == 0}')
    stage(v, '1 congruences matched')
    # handles
    n1 = 6788513 * (v[16742] - v[19083])
    n2 = 14572767 * (v[14681] - v[24483])
    print(f'  26731 numerator / p exact? {n1 % P == 0}')
    print(f'  x_27571 = p ? {v[27571] == P}   9193 numerator / p exact? {n2 % P == 0}')
    seeds = {}
    if n1 % P == 0: seeds[33787] = n1 // P
    if v[27571] == P and n2 % P == 0: seeds[19890] = n2 // P
    ripple(v, seeds)
    nz, f = stage(v, '2 handles set')
    for a in nz:
        print(f'    atom {a} ({len(a2e.get(a,[]))} eqs) gate={atom_out.get(a)}: {src[a][:100]}')
    H.save_assignment(v, 'zero9_out.json')
