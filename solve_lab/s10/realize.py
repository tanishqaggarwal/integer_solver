"""S11 step 15: aim rho = R1/R2 at a 6-subset's required ratio and BUILD the state.

Each 6-subset S has a 1-dimensional kernel w, so A = lambda*w, and the two
congruences become
    lambda*(w2+w3) == R1 ,  lambda*(w5-w4) == R2   (mod p)
solvable for lambda iff  R1*beta == R2*alpha (mod p).
R1 = 5113045*x_7075*x_9118 and R2 = x_8731 are set by two free inputs, so aim
x_9118 at the residue that makes it hold, then realise A exactly and measure.
"""
import os, sys, json, math
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
from frame2 import definer, ORDER, FREE, CHECKS, fwd, score
P = ad.P
SEVEN = [22229, 22230, 35758, 35759, 35760, 35761, 35762]
J = json.load(open(os.path.join(HERE, 'rhs.json')))
E, rows = J['E'], J['rows']
cands = [(tuple(S), w, al, be) for S, w, al, be in J['found']]
base = L.load(os.path.join(LAB, 'best', 'new_instance_partial_39026.json'))
fwd(base)
print(f'delivered witness {score(base)}')

def realize(v, A, x9118):
    v = list(v)
    v[9118] = x9118
    A0, A1, A2, A3, A4, A5, A6 = A
    R1 = 5113045 * v[7075] * v[9118]
    v[29854] = R1 - A3
    if (v[29854] - A2) % P: return None
    v[1329] = (v[29854] - A2) // P
    v[31864] = A5 - v[7075] * v[8731]
    if (v[31864] - A4) % P: return None
    v[10903] = (v[31864] - A4) // P
    v[642] = A6 + v[17325] * P
    v[28730] = A1 + v[9413] * P
    fwd(v, rounds=6)
    v[7068] = A0 + v[2099] + 7376877 * v[642]
    fwd(v, rounds=6)
    v[7068] = A0 + v[2099] + 7376877 * v[642]
    return v

R2 = base[8731] % P
c5113 = 5113045 * base[7075] % P
best = (score(base), None)
tried = 0
seenrho = {}
for S, w, al, be in cands:
    if al % P == 0 or be % P == 0: continue
    rho = (al % P) * pow(be % P, -1, P) % P        # required R1/R2
    R1 = rho * R2 % P
    x9 = R1 * pow(c5113, -1, P) % P                # required x_9118 mod p
    lam = R1 * pow(al % P, -1, P) % P              # lambda mod p
    if (lam * (be % P) - R2) % P: continue
    key = (rho, lam)
    if key in seenrho: continue
    seenrho[key] = 1
    # keep x_9118 near its old magnitude
    k = (base[9118] - x9) // P
    x9118 = x9 + k * P
    A = [lam * c for c in w]
    v = realize(base, A, x9118)
    if v is None: continue
    tried += 1
    av = L.all_atom_values(v)
    fail = L.failing_eqs(av)
    s = L.NEQ - len(fail)
    nz = [a for a in range(L.NA) if av[a]]
    if s > best[0]:
        best = (s, S)
        T.save(v, os.path.join(HERE, f'REAL_{s}.json'))
        print(f'  *** {s} on eqs {[E[i] for i in S]}  nonzero {nz}', flush=True)
    if tried <= 6:
        print(f'  eqs {[E[i] for i in S]}: score {s}  nonzero {len(nz)} '
              f'{nz[:8]}', flush=True)
print(f'\ntried {tried} distinct (rho, lambda); BEST {best[0]} on {best[1]}')
