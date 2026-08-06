"""The all-zero branch:  x_2081 = 0 turns the MUX fully off so x_2099 = x_19964 = 0,
and x_24601 = 0 drives the circuit's output wire x_12186 to 0 as well.
Those two requirements are CONSISTENT (both sides 0), unlike every other quadrant.
Build it and close the pins that the two flips light up."""
import pickle, sys
import harness as H
exec(open('repair.py').read().split('if __name__')[0])
P = 2**256 - 2**32 - 977
d = pickle.load(open('atoms.pkl', 'rb')); src = d['atom_src']
a2e = {a: set(e) for a, e in pickle.load(open('atom2eq.pkl', 'rb')).items()}


def allnz(v):
    return sorted(a for a, Pp in enumerate(polys) if evalpoly(Pp, v) != 0)


def stage(v, tag, show=0):
    nz = allnz(v)
    codes, _ = H.load_equations(); f = H.evaluate(codes, v)
    print(f'[{tag}] atoms={len(nz)}  EQ {len(codes)-len(f)}/{len(codes)} ({len(f)} failing)')
    if show:
        for a in nz[:show]:
            print(f'    atom {a} ({len(a2e.get(a,[]))} eqs) gate={atom_out.get(a)}: {src[a][:95]}')
    return f, nz


if __name__ == '__main__':
    v = H.load_assignment('../best/new_instance_partial_39022.json')
    ripple(v, {2081: 0, 24601: 0})
    stage(v, '1 both control bits 0')
    # x_2081=0 unpins x_6418/x_12553 (atoms 3576/3578 vanish) -> set them to 0
    # x_24601=0 lights x_28180=1 -> atoms 1048/1050 pin x_22152, x_33462 to 0 mod p
    ripple(v, {6418: 0, 12553: 0, 22152: 0, 33462: 0})
    stage(v, '2 unpinned data + complementary pins')
    # the circuit output wire and the verified value both go to 0
    ripple(v, {22649: 0, 14853: 0})
    stage(v, '3 output wire = 0')
    ok, hist = repair_loop(v, rounds=30, verbose=False)
    f, nz = stage(v, '4 after repair', show=25)
    print('  history:', hist)
    print('  core: x_15298=%d u=x_29322=%d w=x_3558=%d' % (v[15298], v[29322], v[3558]))
    print('  atom 22229 (A) =', evalpoly(polys[22229], v))
    print('  atom 22231 (B) =', evalpoly(polys[22231], v))
    H.save_assignment(v, 'zero_out.json')
