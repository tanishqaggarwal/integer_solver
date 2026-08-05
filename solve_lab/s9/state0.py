"""Build S0 = DAG-consistent state (strategy A seeds) and dump residuals."""
import pickle, sys, json
import harness as H
exec(open('repair.py').read().split('if __name__')[0])
P = 2**256-2**32-977
NV = 38748
checks = [a for a in range(len(polys)) if a not in atom_out]
freeinp = [x for x in range(NV) if x not in definer]

def make_S0():
    v = H.load_assignment('../best/new_instance_partial_39022.json')
    ripple(v, {7068: v[2099] + 7376877*v[642], 4432: v[19964] + v[28730]})
    return v

if __name__ == '__main__':
    v = make_S0()
    H.save_assignment(v, 'S0.json')
    nz = [(a, evalpoly(polys[a], v)) for a in checks if evalpoly(polys[a], v) != 0]
    print('S0 nonzero checks:')
    for a, val in nz: print(f'  atom {a}: {val}   %p={val%P}')
    # affineness test on x_24548
    for k in (1, 2, 1000):
        v2 = list(v); ripple(v2, {24548: v[24548]+k})
        d = {a: evalpoly(polys[a], v2)-evalpoly(polys[a], v) for a in checks
             if evalpoly(polys[a], v2) != evalpoly(polys[a], v)}
        print(f'  x_24548 += {k}: ' + ', '.join(f'{a}:{dv//k if dv%k==0 else str(dv)+"(nonlin)"}' for a,dv in sorted(d.items())))
