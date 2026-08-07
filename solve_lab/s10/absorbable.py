"""S11 step 6: which checks only need mod-p preservation?

A check with a solo free handle h (granularity exactly p) absorbs any change that
is a multiple of p.  Those rows need only  delta-residue == 0 (mod p).  Every
OTHER check must hold exactly -- and those are precisely the rows whose response
to large moves is nonlinear, so a linear veto from them is not trustworthy.
Split the closure on that line.
"""
import os, sys, collections, json
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
P = ad.P
definer, atom_out = L.definer, L.atom_out
FREE = set(ad.FREE)
v = L.load(os.path.join(HERE, 'mod9118_0.json'))
CHECKS = [a for a in range(L.NA) if a not in atom_out]
solo = [h for h in FREE if len(L.var_atoms[h]) == 1]
print(f'solo free handles: {len(solo)}')

def dz(a, w):
    """exact integer d(atom a)/d(x_w) at v."""
    s = 0
    for m, c in L.polys[a].items():
        k = m.count(w)
        if k == 0: continue
        if k == 1:
            t = c
            for z in m:
                if z != w: t *= v[z]
            s += t
        else:
            s += 2 * c * v[w]
    return s

absorb = collections.defaultdict(list)     # check -> handles that move it by k*p
gran = collections.Counter()
for h in solo:
    a0 = L.var_atoms[h][0] if isinstance(L.var_atoms[h], list) else list(L.var_atoms[h])[0]
    d = dz(a0, h)
    if d == 0: gran['dormant'] += 1; continue
    gran['p' if d % P == 0 else 'other'] += 1
    if d % P: continue
    # a0 defines some variable t; the change propagates to t's consumers
    if a0 in atom_out:
        t = atom_out[a0][1]
        dt = dz(a0, t)
        if dt == 0: continue
        for c in L.var_atoms[t]:
            if c != a0 and c in set(CHECKS): absorb[c].append(h)
    else:
        absorb[a0].append(h)
print(f'handle granularity: {dict(gran)}')
ABS = set(absorb)
print(f'checks with a p-absorbing solo handle: {len(ABS)} of {len(CHECKS)}')
av = L.all_atom_values(v)
print(f'   the two failing ones absorbable? '
      f'a21617 {21617 in ABS}, a29539 {29539 in ABS}')
json.dump({'absorbable': sorted(ABS)}, open(os.path.join(HERE, 'absorbable.json'), 'w'))
print('saved absorbable.json')
