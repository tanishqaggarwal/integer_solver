"""S10 step 94: the missing move class -- CONGRUENCE moves on free inputs.

a35759 = -x_29854 + 5113045*x_7075*x_9118 with x_29854 = x_1329*p (a35758).
No value of any single variable zeroes it; what it needs is  p | 5113045*x_9118,
i.e. x_9118 == 0 (mod p).  x_9118 is a FREE INPUT, so that is simply a choice --
but no search so far could ever propose it, because rounding x_9118 to a multiple
of p does not zero any atom on its own.
"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
from repair2 import score, candidates
P = ad.P
v0 = L.load(os.path.join(HERE, 'br10.json'))
print(f'br10 score {score(v0)}; x_9118 mod p = {str(v0[9118] % P)[:30]}...')
av = L.all_atom_values(v0)
print(f'a35759 = {str(av[35759])[:40]}   a35758 = {av[35758]}')
print(f'x_1329={str(v0[1329])[:24]} x_22665={str(v0[22665])[:24]} (== p ? {v0[22665]==P})')

for k in (0, 1, -1):
    v = list(v0)
    q = v[9118] // P
    v[9118] = (q + k) * P
    ad.fwd(v, rounds=6)
    a = L.all_atom_values(v)
    nz = [b for b in range(L.NA) if a[b]]
    print(f'\nx_9118 <- p*({q}+{k}) : score {score(v)}  nonzero {nz}')
    # now the second half of the move: close a35759 through x_1329
    for u in (1329,):
        nv = T.solve_lin(35758, u, [x if i != 29854 else 5113045 * v[7075] * v[9118]
                                    for i, x in enumerate(v)])
        print(f'   x_1329 target = {str(nv)[:30] if nv is not None else None}')
        if nv is not None:
            w = list(v); w[u] = nv
            ad.fwd(w, rounds=6)
            aw = L.all_atom_values(w)
            nzw = [b for b in range(L.NA) if aw[b]]
            print(f'   after x_1329 : score {score(w)}  nonzero {nzw}')
            if score(w) >= score(v0):
                T.save(w, os.path.join(HERE, f'mod9118_{k}.json'))
                print(f'   saved mod9118_{k}.json')
