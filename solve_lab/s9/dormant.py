"""Dormant product gates: atoms containing a monomial x_i*x_j where BOTH factors are currently
0 mod p.  Single-variable perturbation says neither factor moves anything (0*d = 0), so every
census and Jacobian in this campaign is blind to them -- but moving BOTH injects d_i*d_j.

This is second-order freedom that no linear analysis can see.  Find it."""
import pickle, sys, collections
import harness as H
exec(open('repair.py').read().split('if __name__')[0])
P = 2**256 - 2**32 - 977
d = pickle.load(open('atoms.pkl', 'rb')); src = d['atom_src']
a2e = {a: set(e) for a, e in pickle.load(open('atom2eq.pkl', 'rb')).items()}
boolv = set(pickle.load(open('boolvars.pkl', 'rb')))
NV = 38748
freeinp = set(x for x in range(NV) if x not in definer)


def dormant_gates(v):
    """(atom, i, j, coeff) for monomials x_i*x_j with both factors 0 mod p"""
    out = []
    for a, Pp in enumerate(polys):
        for m, c in Pp.items():
            if len(m) == 2 and m[0] != m[1]:
                i, j = m
                if v[i] % P == 0 and v[j] % P == 0:
                    out.append((a, i, j, c))
    return out


if __name__ == '__main__':
    src_file = sys.argv[1] if len(sys.argv) > 1 else 'chase_out.json'
    v = H.load_assignment(src_file)
    dg = dormant_gates(v)
    print(f'dormant product monomials (both factors 0 mod p): {len(dg)}')
    # keep those where at least one factor is a FREE input and neither is boolean-locked
    usable = [(a, i, j, c) for a, i, j, c in dg
              if (i in freeinp or j in freeinp) and not (i in boolv and j in boolv)]
    print(f'  with at least one FREE factor and not both boolean: {len(usable)}')
    # of those, the ones whose atom is a GATE (so activating it moves a real wire)
    gates = [(a, i, j, c) for a, i, j, c in usable if a in atom_out]
    print(f'  of which are gate atoms (activating moves a downstream wire): {len(gates)}')
    both_free = [(a, i, j, c) for a, i, j, c in usable if i in freeinp and j in freeinp]
    print(f'  with BOTH factors free inputs: {len(both_free)}')
    for a, i, j, c in both_free[:20]:
        t = atom_out.get(a)
        print(f'    atom {a} out={t} : {c}*x_{i}*x_{j}   {src[a][:70]}')
    pickle.dump({'dormant': dg, 'usable': usable, 'gates': gates, 'both_free': both_free},
                open('dormant.pkl', 'wb'))
