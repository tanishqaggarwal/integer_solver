"""On the x_4287=1 MUX branch the three newly-lit loads (x_2239, x_31731, x_9106) are all
integer combinations of two base quantities x_4306 and x_27177, exactly like L1/L2/L3 are of S,T.
So the branch's obligation is  x_4306 = 0  and  x_27177 = 0  (mod p).
Measure their dependence on the free knobs and solve."""
import pickle, sys, itertools
import harness as H
exec(open('repair.py').read().split('if __name__')[0])
P = 2**256 - 2**32 - 977
roots = pickle.load(open('roots.pkl', 'rb'))
checks = [a for a in range(len(polys)) if a not in atom_out]
rp = {a: (roots[a] if a in roots else polys[a]) for a in checks}

KNOBS = [9118, 8731, 14865, 31861, 6418, 12553]
TARGETS = [4306, 27177, 2239, 31731, 9106]

if __name__ == '__main__':
    v0 = H.load_assignment('mux_out.json')
    print('base residues mod p:')
    for t in TARGETS:
        print(f'   x_{t} % p = {v0[t] % P}')
    # verify the 3-from-2 claim numerically
    print('\ncheck  x_2239 = 3494591*x_27177 + 14240157*x_4306 :',
          (v0[2239] - 3494591*v0[27177] - 14240157*v0[4306]) % P == 0)
    print('check  x_31731 = 15964591*x_27177 + 13881285*x_4306 :',
          (v0[31731] - 15964591*v0[27177] - 13881285*v0[4306]) % P == 0)
    print('check  x_9106  = 7204959*x_27177 + 6822253*x_4306  :',
          (v0[9106] - 7204959*v0[27177] - 6822253*v0[4306]) % P == 0)
    print('\nsensitivity of (x_4306, x_27177) mod p to each free knob (delta per +1):')
    sens = {}
    for k in KNOBS:
        v = list(v0); ripple(v, {k: v0[k] + 1})
        d = ((v[4306]-v0[4306]) % P, (v[27177]-v0[27177]) % P)
        sens[k] = d
        print(f'   x_{k}: d4306={d[0]}  d27177={d[1]}')
    # linearity check
    print('\nlinearity check (knob += 2 gives twice the delta):')
    for k in KNOBS[:3]:
        v = list(v0); ripple(v, {k: v0[k] + 2})
        d = ((v[4306]-v0[4306]) % P, (v[27177]-v0[27177]) % P)
        print(f'   x_{k}: {d == ((2*sens[k][0]) % P, (2*sens[k][1]) % P)}')
    # try to solve  A*z = -(base)  mod p using pairs of knobs
    print('\nsolving x_4306 = x_27177 = 0 mod p with knob pairs:')
    b = ((-v0[4306]) % P, (-v0[27177]) % P)
    for k1, k2 in itertools.combinations(KNOBS, 2):
        a11, a21 = sens[k1]; a12, a22 = sens[k2]
        det = (a11*a22 - a12*a21) % P
        if det == 0: continue
        inv = pow(det, P-2, P)
        z1 = (b[0]*a22 - b[1]*a12) % P * inv % P
        z2 = (a11*b[1] - a21*b[0]) % P * inv % P
        v = list(v0); ripple(v, {k1: v0[k1] + z1, k2: v0[k2] + z2})
        ok = (v[4306] % P == 0 and v[27177] % P == 0)
        if ok:
            print(f'   knobs (x_{k1}, x_{k2}): SOLVED  -> x_4306%p={v[4306]%P} x_27177%p={v[27177]%P}')
            H.save_assignment(v, f'mux2_{k1}_{k2}.json')
            nz = sorted(a for a, Pp in enumerate(polys) if evalpoly(Pp, v) != 0)
            print(f'      nonzero atoms now: {nz}')
