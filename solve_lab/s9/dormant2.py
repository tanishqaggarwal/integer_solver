"""Strictly dormant product gates at the 39,026 witness: monomials x_i*x_j with BOTH factors
EXACTLY 0 (not merely 0 mod p -- the p-wire is 0 mod p but is the ordinary handle mechanism).

Single-variable analysis is blind to these: moving either factor alone gives 0*d = 0.  Moving both
injects d_i*d_j.  Test whether any of them can move the two surviving congruences D1, D2 mod p,
which every linear census says are immovable.
"""
import pickle, sys, time
import harness as H
exec(open('repair.py').read().split('if __name__')[0])
P = 2**256 - 2**32 - 977
d = pickle.load(open('atoms.pkl', 'rb')); src = d['atom_src']
a2e = {a: set(e) for a, e in pickle.load(open('atom2eq.pkl', 'rb')).items()}
boolv = set(pickle.load(open('boolvars.pkl', 'rb')))
NV = 38748
freeinp = set(x for x in range(NV) if x not in definer)
CODES, _ = H.load_equations()


def resid(v):
    return ((v[7068] - v[2099] - 7376877 * v[642]) % P,
            (v[4432] - v[19964] - v[28730]) % P)


def strict_dormant(v):
    seen = set(); out = []
    for a, Pp in enumerate(polys):
        for m, c in Pp.items():
            if len(m) == 2 and m[0] != m[1]:
                i, j = m
                if v[i] == 0 and v[j] == 0 and (i, j) not in seen:
                    seen.add((i, j)); out.append((a, i, j, c))
    return out


if __name__ == '__main__':
    v0 = H.load_assignment('../best/new_instance_partial_39026.json')
    b0 = resid(v0)
    print(f'39,026 witness: D1%p={b0[0]}\n                D2%p={b0[1]}')
    dg = strict_dormant(v0)
    print(f'\nstrictly dormant monomials (both factors exactly 0): {len(dg)}')
    # keep pairs where both factors are settable: free inputs, and not both boolean
    cand = []
    seenp = set()
    for a, i, j, c in dg:
        if (i, j) in seenp: continue
        if i in freeinp and j in freeinp:
            seenp.add((i, j)); cand.append((a, i, j, c))
    print(f'  with BOTH factors free inputs: {len(cand)}')
    print('\nactivating each and measuring the effect on (D1, D2) mod p ...')
    hits = []
    t0 = time.time()
    for n, (a, i, j, c) in enumerate(cand):
        w = list(v0)
        try:
            ripple(w, {i: 1, j: 1})
        except Exception:
            continue
        b = resid(w)
        if b != b0:
            hits.append((a, i, j, (b[0] - b0[0]) % P, (b[1] - b0[1]) % P))
        if n % 200 == 0:
            print(f'   ...{n}/{len(cand)}  {time.time()-t0:.0f}s  hits={len(hits)}', flush=True)
    print(f'\ndormant gates whose activation MOVES (D1,D2) mod p: {len(hits)}')
    for a, i, j, d1, d2 in hits[:20]:
        print(f'   atom {a}: x_{i}*x_{j}  dD1={str(d1)[:26]}  dD2={str(d2)[:26]}')
    pickle.dump(hits, open('dormant_hits.pkl', 'wb'))
