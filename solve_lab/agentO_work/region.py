"""The 39,026 residual region has SEVEN variables private to it (they occur in no atom
   outside the 8-atom region): x_642, x_1329, x_9413, x_10903, x_17325, x_29854, x_31864.
   The 12 equations that touch the region have rank 7 over the 8 residuals, with kernel
   spanned by (a36662 = +t, a36663 = -t) because those two atoms have IDENTICAL coefficient
   columns in all 12.  So a full solve of the region needs only:
       a23616 = a23617 = a36659 = a36660 = a36661 = a36664 = 0  and  a36663 = -a36662.
   This script computes the required private values exactly and tests divisibility.
"""
import sys, json, collections, itertools
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentO_work')
import simO, engine as E, harness as H

OD = '/home/user/integer_solver/solve_lab/agentO_work'
d = json.load(open('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json'))
v = [0] * E.NV
for k, x in d.items():
    v[int(k.split('_')[1])] = int(x)

x7068, x2099, x28599, x28730, x17499 = v[7068], v[2099], v[28599], v[28730], v[17499]
x22665, x7075, x9118, x28961, x8731 = v[22665], v[7075], v[9118], v[28961], v[8731]
print('x_7075 =', x7075)
J = x7075 * x8731                      # a36662, fixed
print('a36662 = x_7075*x_8731, bits =', abs(J).bit_length())

need = {}
ok = True


def dv(name, num, den):
    global ok
    if den == 0 or num % den:
        print(f'  {name}: NOT divisible  ({abs(num).bit_length()}b / {abs(den).bit_length()}b)')
        ok = False
        return None
    q = num // den
    print(f'  {name} = {abs(q).bit_length()}b  OK')
    return q


print('required private values for a FULL region solve:')
x642 = dv('x_642 from a23616 : (x_7068-x_2099)/7376877', x7068 - x2099, 7376877)
if x642 is not None:
    need[642] = x642
    x17325 = dv('x_17325 from a36664: x_642/x_28599', x642, x28599)
    if x17325 is not None:
        need[17325] = x17325
x9413 = dv('x_9413 from a23617 : x_28730/x_17499', x28730, x17499)
if x9413 is not None:
    need[9413] = x9413
x29854 = 5113045 * (x7075 * x9118)
need[29854] = x29854
print(f'  x_29854 from a36660: 5113045*x_7075*x_9118 = {abs(x29854).bit_length()}b  OK')
x1329 = dv('x_1329 from a36659 : x_29854/x_22665', x29854, x22665)
if x1329 is not None:
    need[1329] = x1329
x31864 = -J
need[31864] = x31864
print(f'  x_31864 from a36663: -a36662 = {abs(x31864).bit_length()}b  OK')
x10903 = dv('x_10903 from a36661: x_31864/x_28961', x31864, x28961)
if x10903 is not None:
    need[10903] = x10903

print('\nall divisibilities hold:', ok)
w = list(v)
for u, val in need.items():
    w[u] = val
bad = E.badatoms(w)
ff = E.eqfails(bad)
print('after substitution: bad atoms =', sorted(bad), ' failing eqs =', len(ff), sorted(ff))
print('score =', 39033 - len(ff))
json.dump({f"x_{i}": str(int(w[i])) for i in range(E.NV) if w[i] != 0},
          open(f'{OD}/region_try_{39033-len(ff)}.json', 'w'))
print('wrote', f'{OD}/region_try_{39033-len(ff)}.json')
