"""Full construction on the x_4287 = 1 MUX branch:
  1. flip the MUX control so x_2099/x_19964 come from the FREE inputs x_9118/x_8731
  2. close the pins the flip lights (x_31861, x_14865)
  3. solve the branch obligation  x_4306 = x_27177 = 0 (mod p)  with (x_9118, x_8731)
  4. set the three quotient handles so the newly-lit loads vanish exactly in Z
  5. repair the mirrors, keeping the core control difference u = 0
"""
import pickle, sys, itertools
import harness as H
exec(open('repair.py').read().split('if __name__')[0])
P = 2**256 - 2**32 - 977
d = pickle.load(open('atoms.pkl', 'rb')); src = d['atom_src']
a2e = {a: set(e) for a, e in pickle.load(open('atom2eq.pkl', 'rb')).items()}
H1 = 119562606790549640390870952418684367882170154220603339634805704742270834564330392414192110
H2 = 113141528427610260107049117992526537105383080782811760722361109500341947028737388716982706


def allnz(v):
    return sorted(a for a, Pp in enumerate(polys) if evalpoly(Pp, v) != 0)


def stage(v, tag, show=False):
    nz = allnz(v)
    codes, _ = H.load_equations(); f = H.evaluate(codes, v)
    print(f'[{tag}] atoms={len(nz)} {nz}   EQ {len(codes)-len(f)}/{len(codes)} ({len(f)} failing)')
    if show:
        for a in nz:
            print(f'    atom {a} ({len(a2e.get(a,[]))} eqs) gate={atom_out.get(a)}: {src[a][:95]}')
    return f, nz


def build():
    v = H.load_assignment('../best/new_instance_partial_39022.json')
    ripple(v, {4287: 1})
    ripple(v, {9118: v[7068] - 7376877*v[642], 8731: v[4432] - v[28730]})
    ripple(v, {31861: H1, 14865: H2})
    stage(v, '1 mux+pins')
    # --- 3. solve the branch obligation with (x_9118, x_8731)
    v0 = list(v)
    sens = {}
    for k in (9118, 8731):
        w = list(v0); ripple(w, {k: v0[k] + 1})
        sens[k] = ((w[4306]-v0[4306]) % P, (w[27177]-v0[27177]) % P)
    b = ((-v0[4306]) % P, (-v0[27177]) % P)
    a11, a21 = sens[9118]; a12, a22 = sens[8731]
    det = (a11*a22 - a12*a21) % P
    inv = pow(det, P-2, P)
    z1 = (b[0]*a22 - b[1]*a12) % P * inv % P
    z2 = (a11*b[1] - a21*b[0]) % P * inv % P
    ripple(v, {9118: v0[9118] + z1, 8731: v0[8731] + z2})
    print(f'   x_4306%p={v[4306]%P}  x_27177%p={v[27177]%P}')
    stage(v, '2 obligation solved')
    # --- 4. set the quotient handles for the newly-lit loads
    seeds = {}
    n = 6122989*v[2239]
    if n % P == 0: seeds[6947] = n // P
    n = -v[31731]
    if n % P == 0: seeds[33168] = n // P
    n = v[9106]
    if n % (13523997*P) == 0: seeds[950] = n // (13523997*P)
    print('   handle seeds:', {k: str(x)[:14]+'...' for k, x in seeds.items()})
    ripple(v, seeds)
    stage(v, '3 handles set')
    # --- 5. mirrors: x_14853 = x_1308 and x_24548 = x_25442 (both free inputs)
    ripple(v, {14853: v[1308], 24548: v[25442]})
    stage(v, '4 mirrors closed', show=True)
    print('   core: x_15298=%d  u=x_29322=%d  w=x_3558=%d' % (v[15298], v[29322], v[3558]))
    ok, hist = repair_loop(v, rounds=15, verbose=False)
    stage(v, '5 after repair', show=True)
    return v


if __name__ == '__main__':
    v = build()
    H.save_assignment(v, 'mux3_out.json')
