"""S10 step 59: price the forced OR gate.

Every load pin  bit*(x_B - HUGE) - s*x_C  is satisfied FOR FREE when its bit is 0
(the pin collapses to -s*x_C and the handle takes 0).  The HUGE constants -- the
only source of values that are not multiples of p -- therefore enter ONLY through
set bits.  Session 9 established that the OR gate x_9274 = 1 is forced and can be
satisfied only by a gate-bit.  But that is an ATOM-space statement.  In equation
space the question is simply: what does violating it COST?

Budget to beat: 7 equations.
"""
import os, sys, collections, json
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L
import tools as T
import ad
from newton import BOOL

P = ad.P
atom_out = L.atom_out
base = L.load(os.path.join(HERE, 'forward_state.json'))

print('=== x_9274 and everything that constrains it ===')
d = L.definer.get(9274)
print(f'x_9274 = {base[9274]}   definer: {("a%d" % d) if d is not None else "FREE"}')
if d is not None:
    print(f'   a{d}: {L.atom_src[d][:200]}   price={len(L.atom2eq.get(d,{}))}')
for a in sorted(L.var_atoms[9274]):
    out = L.atom_out.get(a)
    print(f'   a{a:<7} {"GATE->x_%d" % out[1] if out else "CHECK":<16} '
          f'price={len(L.atom2eq.get(a,{})):<4} {L.atom_src[a][:130]}')

print('\n=== all bare pins  (atom = single variable - constant)  and their prices ===')
bare = []
for a in range(L.NA):
    poly = L.polys[a]
    if len(poly) == 2:
        lin = [m for m in poly if len(m) == 1]
        con = [m for m in poly if len(m) == 0]
        if len(lin) == 1 and len(con) == 1:
            bare.append((len(L.atom2eq.get(a, {})), a, lin[0][0], poly[con[0]]))
    elif len(poly) == 1:
        lin = [m for m in poly if len(m) == 1]
        if len(lin) == 1:
            bare.append((len(L.atom2eq.get(a, {})), a, lin[0][0], 0))
bare.sort()
print(f'bare single-variable pins: {len(bare)}')
hist = collections.Counter(p for p, a, u, c in bare)
print('  price histogram:', dict(sorted(hist.items())[:10]))
print('  cheapest:')
for p, a, u, c in bare[:14]:
    tag = 'BOOL' if u in BOOL else ''
    print(f'    a{a:<7} price={p:<4} x_{u:<7} = {str(-c)[:34]:<36} {tag}')

print('\n=== the boolean free inputs currently set to 1 ===')
ones = [u for u in ad.FREE if u in BOOL and base[u] == 1]
print(f'boolean free inputs = 1: {len(ones)} -> {ones[:40]}')

print('\n=== experiment: set every boolean free input to 0 ===')
v = list(base)
for u in ad.FREE:
    if u in BOOL:
        v[u] = 0
ad.fwd(v, rounds=4)
av = L.all_atom_values(v)
nz = [a for a in range(L.NA) if av[a]]
nzc = [a for a in nz if a not in atom_out]
fail = L.failing_eqs(av)
print(f'  nonzero atoms={len(nz)} (checks {len(nzc)}) failing={len(fail)} '
      f'score={L.NEQ-len(fail)}')
print(f'  x_9274 = {v[9274]}')
tot = 0
for a in nzc[:25]:
    pr = len(L.atom2eq.get(a, {}))
    tot += pr
    print(f'    a{a:<7} price={pr:<4} {L.atom_src[a][:105]}')
print(f'  ... {len(nzc)} checks, price sum of shown {tot}')
T.save(v, os.path.join(HERE, 'allbits0.json'))
