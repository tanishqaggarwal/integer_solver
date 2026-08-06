"""Audit every link in the chain that pins x_7068 to CONST.

Chain:  pin 31670 -> x_22152 -> x_29524 -> (atom 2423) -> x_22649 -> x_12186
        -> (core u=0) -> x_14853 -> (mirror 29539) -> x_1308 = x_7068 -> atom 22229 = A

Each link is either an identity/sum gate (rigid) or a congruence check k*(X-Y) = handle.
A congruence link is only rigid if its handle is a multiple of the p-wire.  Find any link whose
handle is NOT p-quantised -- that link can absorb an arbitrary shift and breaks the chain.
"""
import pickle, sys
import harness as H
exec(open('repair.py').read().split('if __name__')[0])
P = 2**256 - 2**32 - 977
d = pickle.load(open('atoms.pkl', 'rb')); src = d['atom_src']
occ = pickle.load(open('atom_occ_all.pkl', 'rb'))
a2e = {a: set(e) for a, e in pickle.load(open('atom2eq.pkl', 'rb')).items()}
NV = 38748
freeinp = set(x for x in range(NV) if x not in definer)

LINKS = [31670, 2423, 22772, 26729, 18869, 29539, 22229, 1048]


def handle_granularity(v, Hv):
    """what values can the handle variable Hv take?  returns a description + the modulus"""
    hd = definer.get(Hv)
    if hd is None:
        return 'FREE (any integer)', 1
    Pp = polys[hd]
    prods = [(m, c) for m, c in Pp.items() if len(m) == 2 and Hv not in m]
    if not prods:
        return f'defined by atom {hd} (non-product): {src[hd][:70]}', None
    m, c = prods[0]
    w1, w2 = m
    g = []
    for a_, b_ in ((w1, w2), (w2, w1)):
        if v[a_] == P: g.append((f'x_{a_} = p (wire)', P))
        elif v[a_] == 0: g.append((f'x_{a_} = 0', 0))
        else: g.append((f'x_{a_} = {str(v[a_])[:14]}...', v[a_]))
    free_side = [x for x in (w1, w2) if x in freeinp]
    mods = [v[w1], v[w2]]
    gran = None
    if v[w1] == P and w2 in freeinp: gran = P
    elif v[w2] == P and w1 in freeinp: gran = P
    return (f'{Hv} = x_{w1}*x_{w2}  [{g[0][0]} , {g[1][0]}]  free={free_side}'), gran


if __name__ == '__main__':
    v = H.load_assignment('../best/new_instance_partial_39026.json')
    print('link audit at the 39,026 witness\n' + '=' * 70)
    for a in LINKS:
        Pp = polys[a]
        lin = {m[0]: c for m, c in Pp.items() if len(m) == 1}
        val = evalpoly(Pp, v)
        print(f'\natom {a}  ({len(a2e.get(a,[]))} eqs)  value={val}')
        print(f'   {src[a][:150]}')
        # candidate handle = a linear term whose variable is defined by a product
        for t, c in lin.items():
            hd = definer.get(t)
            if hd is None: continue
            hp = polys[hd]
            if any(len(m) == 2 and t not in m for m in hp):
                desc, gran = handle_granularity(v, t)
                print(f'   handle candidate x_{t} (coef {c}): {desc}')
                print(f'      -> granularity {"p" if gran == P else gran}')
