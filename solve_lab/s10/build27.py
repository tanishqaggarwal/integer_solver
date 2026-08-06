"""S10 step 69: realise the 8-atom target and MEASURE it with the real checker.

Target A = (a22229,a22230,a35758,a35759,a35760,a35761,a35762,a22231):
    x_642    = (D - A1)/7376877        x_17325 = (x_642 - A7)/p
    x_28730  = K - A8                  x_9413  = (x_28730 - A2)/p
    x_1329 chosen so 5113045 | (A3+A4+p*x_1329); x_29854 = A3 + p*x_1329;
    x_9118   = (A3+A4+p*x_1329)/5113045
    x_10903  = 0, x_31864 = A5, x_8731 = A6 - A5
x_4432 is NOT touched, so nothing downstream moves.
"""
import os, sys, json, math
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L
import tools as T

P = 2**256 - 2**32 - 977
BEST = os.path.join(LAB, 'best', 'new_instance_partial_39026.json')
base = L.load(BEST)
spec = json.load(open(os.path.join(HERE, 'a22231_best.json')))
A = [int(x) for x in spec['A']]
NZ = spec['NZ']
print('target atoms', NZ)
print('target values', [str(x)[:24] for x in A])

BLOCK = set(NZ) | {37887}
D = base[7068] - base[2099]
K = base[4432] - base[19964]
A1, A2, A3, A4, A5, A6, A7, A8 = A

v = list(base)
# free k*p shift of x_7068 (s10/repairD.py) fixes D mod 7376877
k = ((A1 - D) * pow(P, -1, 7376877)) % 7376877
print(f'\nshifting x_7068 by k*p with k = {k}')
L.ripple(v, {7068: base[7068] + k * P}, block=BLOCK)
for a in (29539, 40826):
    opts = []
    for u in sorted(L.avars[a]):
        nv = T.solve_lin(a, u, v)
        if nv is not None and nv != v[u]:
            opts.append((len(L.var_eqs[u]), u, nv))
    opts.sort()
    if opts:
        _, u, nv = opts[0]
        L.ripple(v, {u: nv}, block=BLOCK)
        print(f'   closed a{a} via x_{u}')
D = v[7068] - v[2099]
K = v[4432] - v[19964]
assert (D - A1) % 7376877 == 0, 'x_642 not integral'
x642 = (D - A1) // 7376877
assert (x642 - A7) % P == 0, 'x_17325 not integral'
x17325 = (x642 - A7) // P
x28730 = K - A8
assert (x28730 - A2) % P == 0, 'x_9413 not integral'
x9413 = (x28730 - A2) // P
t = (-(A3 + A4) * pow(P, -1, 5113045)) % 5113045
num = A3 + A4 + P * t
assert num % 5113045 == 0
x9118 = num // 5113045
x29854 = A3 + P * t

seeds = {642: x642, 17325: x17325, 28730: x28730, 9413: x9413,
         1329: t, 29854: x29854, 9118: x9118,
         10903: 0, 31864: A5, 8731: A6 - A5}
print('\nseeds (x_4432 deliberately untouched):')
for k, val in seeds.items():
    print(f'   x_{k:<7} -> {str(val)[:34]}')
ch, _ = L.ripple(v, seeds, block=BLOCK)
print(f'ripple changed {len(ch)} variables; x_4432 moved: {v[4432] != base[4432]}')

av = L.all_atom_values(v)
nz = [a for a in range(L.NA) if av[a]]
fail = L.failing_eqs(av)
print(f'\nrealised atoms: {[(a, str(av[a])[:20]) for a in nz]}')
print(f'target matched: {[av[NZ[i]] == A[i] for i in range(len(NZ))]}')
print(f'nonzero atoms {nz}')
print(f'FAILING {len(fail)}  SCORE {L.NEQ - len(fail)}   (previous best 39026)')
print(f'failing equations: {fail}')
if len(fail) < 7:
    out = os.path.join(LAB, 'best', f'new_instance_partial_{L.NEQ-len(fail)}.json')
    T.save(v, out)
    print(f'\n*** IMPROVEMENT -- saved {out}')
else:
    T.save(v, os.path.join(HERE, f'cand27_{L.NEQ-len(fail)}.json'))
