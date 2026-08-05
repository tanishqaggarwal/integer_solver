"""Close the all-zero branch.
Remaining primitives after zero.py:
  30976  x_15029 = 2104445*x_22162,  x_15029 = x_32762*x_8386 with x_32762 = p  ->  need p | x_22162
  30978  x_36202 = -x_30213,         x_36202 = x_35409*x_21868 with x_35409 = p ->  need p | x_30213
  23000  x_9274 = OR(x_7715, x_34554) and x_9274 = x_2300 = 1 (pinned)          ->  need exactly one bit set
  688    x_18956 = C2const (mod p)
  1618   x_24468 = C1const (mod p)
Both x_22162 and x_30213 are FREE inputs, so the first two close by setting them to 0.
"""
import pickle, sys
import harness as H
exec(open('repair.py').read().split('if __name__')[0])
P = 2**256 - 2**32 - 977
d = pickle.load(open('atoms.pkl', 'rb')); src = d['atom_src']
a2e = {a: set(e) for a, e in pickle.load(open('atom2eq.pkl', 'rb')).items()}
boolv = set(pickle.load(open('boolvars.pkl', 'rb')))
NV = 38748
freeinp = [x for x in range(NV) if x not in definer]
bfree = [b for b in freeinp if b in boolv]


def allnz(v):
    return sorted(a for a, Pp in enumerate(polys) if evalpoly(Pp, v) != 0)


def stage(v, tag, show=0):
    nz = allnz(v)
    codes, _ = H.load_equations(); f = H.evaluate(codes, v)
    print(f'[{tag}] atoms={len(nz)} {nz}  EQ {len(codes)-len(f)}/{len(codes)} ({len(f)} failing)')
    return f, nz


if __name__ == '__main__':
    v = H.load_assignment('zero_out.json')
    stage(v, '0 start')
    ripple(v, {22162: 0, 30213: 0})
    stage(v, '1 zero the two free carriers')
    # find a boolean free input that sets x_7715 = 1 while keeping x_34554 = 0 (so x_15298 stays 0)
    cands = []
    for b in bfree:
        w = list(v); ripple(w, {b: 1 - v[b]})
        if w[9274] == 1 and w[15298] == 0:
            cands.append((b, w[7715], w[34554]))
    print(f'bits giving x_9274 = 1 with the core still off: {len(cands)} -> {[c[0] for c in cands[:12]]}')
    best = None
    for b, a, c in cands:
        w = list(v); ripple(w, {b: 1 - v[b]})
        ok, hist = repair_loop(w, rounds=12, verbose=False)
        nz = allnz(w)
        if best is None or len(nz) < len(best[1]): best = (b, nz, w)
    if best:
        b, nz, w = best
        print(f'best bit x_{b}: atoms={nz}')
        f, nz = stage(w, f'2 with bit x_{b}')
        for a in nz:
            print(f'    atom {a} ({len(a2e.get(a,[]))} eqs) gate={atom_out.get(a)}: {src[a][:95]}')
        H.save_assignment(w, 'zero2_out.json')
