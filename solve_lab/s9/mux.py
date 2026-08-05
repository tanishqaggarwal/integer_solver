"""Flip the MUX control x_4287 -> 1 so that x_2099 and x_19964 are selected from the FREE
inputs x_9118 / x_8731 instead of the pinned constants x_6418 / x_12553, then set them to
satisfy C1 and C2 exactly."""
import pickle, sys
import harness as H
exec(open('repair.py').read().split('if __name__')[0])
P = 2**256 - 2**32 - 977
roots = pickle.load(open('roots.pkl', 'rb'))
checks = [a for a in range(len(polys)) if a not in atom_out]
rp = {a: (roots[a] if a in roots else polys[a]) for a in checks}
d = pickle.load(open('atoms.pkl', 'rb')); src = d['atom_src']
a2e = {a: set(e) for a, e in pickle.load(open('atom2eq.pkl', 'rb')).items()}


def allnz(v):
    return sorted(a for a, Pp in enumerate(polys) if evalpoly(Pp, v) != 0)


def report(v, tag):
    nz = allnz(v)
    codes, _ = H.load_equations(); f = H.evaluate(codes, v)
    print(f'[{tag}] nonzero atoms={nz}')
    print(f'[{tag}] EQUATIONS {len(codes)-len(f)}/{len(codes)}  ({len(f)} failing)')
    return f, nz


if __name__ == '__main__':
    base = sys.argv[1] if len(sys.argv) > 1 else '../best/new_instance_partial_39022.json'
    v = H.load_assignment(base)
    print('before: x_4287=%d x_21279=%d x_2099%%p=%d x_19964%%p=%d'
          % (v[4287], v[21279], v[2099] % P, v[19964] % P))
    ripple(v, {4287: 1})
    print('after flip: x_9062=%d x_21279=%d x_6788=%d x_31033=%d' % (v[9062], v[21279], v[6788], v[31033]))
    print('  x_2099 now = x_9118? ', v[2099] == v[9118], ' x_19964 = x_8731? ', v[19964] == v[8731])
    # C1: x_7068 = x_2099 + 7376877*x_642 ;  C2: x_4432 = x_19964 + x_28730
    want2099 = v[7068] - 7376877 * v[642]
    want19964 = v[4432] - v[28730]
    ripple(v, {9118: want2099, 8731: want19964})
    print('atom 22229 (A) =', evalpoly(polys[22229], v))
    print('atom 22231 (B) =', evalpoly(polys[22231], v))
    f, nz = report(v, 'mux')
    # close the pins that x_4287=1 lights up, then repair the rest
    H1 = 119562606790549640390870952418684367882170154220603339634805704742270834564330392414192110
    H2 = 113141528427610260107049117992526537105383080782811760722361109500341947028737388716982706
    ripple(v, {31861: H1, 14865: H2})
    print('atom 3568 =', evalpoly(polys[3568], v), ' atom 3570 =', evalpoly(polys[3570], v))
    f, nz = report(v, 'mux+pins')
    ok, hist = repair_loop(v, rounds=25, verbose=True)
    f, nz = report(v, 'mux+pins+repair')
    for a in nz:
        print(f'   atom {a} ({len(a2e.get(a,[]))} eqs) gate={atom_out.get(a)}: {src[a][:100]}')
    for a in nz:
        print(f'   atom {a} ({len(a2e.get(a,[]))} eqs) gate={atom_out.get(a)}: {src[a][:110]}')
    H.save_assignment(v, 'mux_out.json')
