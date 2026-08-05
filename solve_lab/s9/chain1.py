"""Kill chain 1: drive the circuit's computed value x_12186 to the pinned constant K1 via the
free input x_22649, moving x_14853 with it so the core's control difference u stays 0."""
import pickle, sys
import harness as H
exec(open('repair.py').read().split('if __name__')[0])
P = 2**256 - 2**32 - 977
roots = pickle.load(open('roots.pkl', 'rb'))
checks = [a for a in range(len(polys)) if a not in atom_out]
rp = {a: (roots[a] if a in roots else polys[a]) for a in checks}
d = pickle.load(open('atoms.pkl', 'rb')); src = d['atom_src']


def nz(v):
    return sorted(a for a, Pp in rp.items() if evalpoly(Pp, v) != 0)


def build(base='../best/new_instance_partial_39022.json', also_29524=True):
    v = H.load_assignment(base)
    K1 = v[2099]                      # the pinned constant (x_2099 = x_6418 = K1)
    delta = K1 - v[12186]
    print(f'delta = K1 - x_12186 = {delta}')
    # x_22152 drives x_29524 1:1, so move it instead of breaking gate 22772
    seeds = {22649: v[22649] + delta, 22152: v[22152] + delta,
             14853: v[14853] + delta, 7068: K1}
    ripple(v, seeds)
    print('after move: x_12186 =', v[12186] == K1, ' x_14853 =', v[14853] == K1,
          ' x_7068 =', v[7068] == K1, ' u = x_29322 =', v[29322])
    print('atom 22229 (A) =', evalpoly(polys[22229], v))
    print('atom 22231 (B) =', evalpoly(polys[22231], v))
    return v


if __name__ == '__main__':
    v = build()
    print('nonzero checks:', nz(v))
    codes, _ = H.load_equations(); f = H.evaluate(codes, v)
    print(f'EQUATIONS: {len(codes)-len(f)}/{len(codes)}  ({len(f)} failing) {f[:20]}')
    H.save_assignment(v, 'chain1_out.json')
