"""S10 step 13: realise a target atom-vector A as an actual variable assignment.

Given A = (a22229,a22230,a35758,a35759,a35760,a35761,a35762) in the achievable set:
  x_642    = (D - A1)/7376877
  x_17325  = (x_642 - A7)/p
  x_28730  = A2,        x_9413  = 0            (then ripple repairs x_4432)
  x_1329  chosen so 5113045 | (A3 + A4 + p*x_1329);  x_29854 = A3 + p*x_1329;
  x_9118   = (A4 + A3 + p*x_1329)/5113045
  x_10903  = 0, x_31864 = A5, x_8731 = A6 - A5
"""
import os, sys, json, math, itertools
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L
import tools as T

P = 2**256 - 2**32 - 977
MOD = 7376877 * P
NZ = [22229, 22230, 35758, 35759, 35760, 35761, 35762]
BEST = os.path.join(LAB, 'best', 'new_instance_partial_39026.json')
BLOCK = set(NZ) | {22231}          # we set these atoms' variables by hand


def build(A, base):
    v = list(base)
    D = v[7068] - v[2099]
    A1, A2, A3, A4, A5, A6, A7 = A
    assert (D - A1) % 7376877 == 0, 'x_642 not integral'
    x642 = (D - A1) // 7376877
    assert (x642 - A7) % P == 0, 'x_17325 not integral'
    x17325 = (x642 - A7) // P
    t = (-(A3 + A4) * pow(P, -1, 5113045)) % 5113045
    num = A3 + A4 + P * t
    assert num % 5113045 == 0
    x9118 = num // 5113045
    x29854 = A3 + P * t
    seeds = {642: x642, 17325: x17325, 28730: A2, 9413: 0,
             1329: t, 29854: x29854, 9118: x9118,
             10903: 0, 31864: A5, 8731: A6 - A5}
    # x_4432 must follow x_28730 so that a22231 stays 0
    seeds[4432] = v[19964] + A2
    L.ripple(v, seeds, block=BLOCK)
    return v


def score(v):
    av = L.all_atom_values(v)
    fail = L.failing_eqs(av)
    return L.NEQ - len(fail), [a for a in range(L.NA) if av[a]], fail


base = L.load(BEST)
print('base score', score(base)[0])

spec = json.load(open(os.path.join(HERE, 'lattice_best.json')))
A = [int(x) for x in spec['A']]
print('target S =', spec['S'])
v = build(A, base)
sc, nz, fail = score(v)
print(f'constructed: score={sc} nz_atoms={nz} failing={fail}')
av = L.all_atom_values(v)
print('realised A =', [av[a] for a in NZ])
print('target   A =', A)
print('match      =', [av[a] for a in NZ] == A)
if sc >= 39027:
    T.save(v, os.path.join(HERE, f'cand_{sc}.json'))
    print(f'saved cand_{sc}.json')
