"""S10 step 4: drive the residual to its minimal form.

From the 7-atom analysis:
  a35762 = x_642   - p*x_17325      -> zero by x_642=0, x_17325=0
  a22230 = x_28730 - p*x_9413       -> zero by x_28730=0, x_9413=0
  a35758 = x_29854 - p*x_1329       -> zero by x_29854=0, x_1329=0
  a35759 = 5113045*x_7075*x_9118 - x_29854 -> zero by x_9118=0
  a35760 = x_31864 - p*x_10903      -> zero by x_31864=0, x_10903=0
  a35761 = x_7075*x_8731 + x_31864  -> zero by x_8731=0
  a22229 = x_7068 - x_2099 - 7376877*x_642  ->  x_7068 == x_2099

So EVERYTHING reduces to the single scalar identity  x_7068 == x_2099.
Build that state, ripple, and measure.
"""
import os, sys, json
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
import lib as L

P = 2**256 - 2**32 - 977
v = L.load(os.path.join(LAB, 'best', 'new_instance_partial_39026.json'))
base_av = L.all_atom_values(v)
print('start: nonzero atoms', [a for a in range(L.NA) if base_av[a]],
      ' failing', len(L.failing_eqs(base_av)))

seeds = {642: 0, 17325: 0, 28730: 0, 9413: 0, 29854: 0, 1329: 0,
         31864: 0, 10903: 0, 9118: 0, 8731: 0}
# block the seven residual atoms from being "repaired" by the ripple: we are
# setting their variables by hand.
BLOCK = {22229, 22230, 35758, 35759, 35760, 35761, 35762}
changed, steps = L.ripple(v, seeds, block=BLOCK)
av = L.all_atom_values(v)
nz = [a for a in range(L.NA) if av[a]]
fail = L.failing_eqs(av)
print(f'after ripple: changed {len(changed)} vars in {steps} steps')
print(f'  nonzero atoms: {len(nz)} {nz}')
for a in nz:
    print(f'    a{a} = {L.atom_src[a][:110]}   -> {str(av[a])[:60]}')
print(f'  failing equations: {len(fail)} -> score {L.NEQ-len(fail)}/{L.NEQ}')
print(f'  x_7068 = {v[7068]}')
print(f'  x_2099 = {v[2099]}')
print(f'  difference = {v[7068]-v[2099]}')
d = v[7068] - v[2099]
print(f'  d mod p        = {d % P}')
print(f'  d mod 7376877  = {d % 7376877}')
json.dump({f'x_{i}': v[i] for i in range(L.NVARS) if v[i] != 0},
          open(os.path.join(HERE, 'clean_state.json'), 'w'))
print('saved clean_state.json')
