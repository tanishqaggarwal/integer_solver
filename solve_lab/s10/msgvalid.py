"""S10 step 86: where IS the linear bit model valid?

msgverify.py showed the sweep's best state is bogus: its two bits were x_2081 and
x_4287, the structural MUX controls, where b*(X - HUGE) has X itself depending on
b so linearity fails.  But the two LARGE groups (75 and 50 bits) are ordinary load
bits whose contribution really is linear.  Validate that empirically, then re-run
the sweep restricted to the groups where the model holds.
"""
import os, sys, json, collections
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L
import tools as T
import ad

P = ad.P
atom_out = L.atom_out
spec = json.load(open(os.path.join(HERE, 'msgsweep.json')))
groups = spec['groups']; CHECKS = spec['checks']
v = L.load(os.path.join(HERE, 'forward_state.json'))
av0 = L.all_atom_values(v)
f0 = L.failing_eqs(av0)
print(f'base failing {len(f0)}')

BIG = [(list(map(int, s)), us) for s, us in groups if len(us) > 1]
SINGLE = [(list(map(int, s)), us) for s, us in groups if len(us) == 1]
print(f'large groups: {[len(us) for s, us in BIG]}; '
      f'singletons: {[us[0] for s, us in SINGLE]}')

print('\n=== validate linearity on the large groups ===')
for s, us in BIG:
    u = us[0]
    w = list(v); w[u] = 1 - w[u]
    ad.fwd(w, rounds=3)
    aw = L.all_atom_values(w)
    ok = True
    for ci, c in enumerate(CHECKS):
        pred = (av0[c] + (1 if v[u] == 0 else -1) * s[ci]) % P
        got = aw[c] % P
        if pred != got:
            ok = False
            print(f'   group of {len(us)}: bit x_{u}, check a{c}: '
                  f'MISMATCH (model {str(pred)[:18]} vs real {str(got)[:18]})')
    nz = [a for a in range(L.NA) if aw[a] and a not in atom_out]
    print(f'   group of {len(us):>3}: bit x_{u} -> model {"MATCHES" if ok else "FAILS"}; '
          f'failing checks {nz}; failing eqs {len(L.failing_eqs(aw))}')

print('\n=== sweep restricted to the two large (linear) groups ===')
sizes = [len(us) for s, us in BIG]
sigs = [s for s, us in BIG]
best = None
hist = collections.Counter()
for c1 in range(sizes[0] + 1):
    for c2 in range(sizes[1] + 1):
        z = 0
        for ci, c in enumerate(CHECKS):
            d = c1 * sigs[0][ci] + c2 * sigs[1][ci]
            if (av0[c] + d) % P == 0:
                z += 1
        hist[z] += 1
        if best is None or z > best[0]:
            best = (z, c1, c2)
print(f'  states swept: {(sizes[0]+1)*(sizes[1]+1)}')
print(f'  histogram of checks zeroed: {dict(sorted(hist.items()))}')
print(f'  BEST: {best[0]} checks zeroed at counts {best[1:]} ')
print('\n=> the 125 ordinary load bits, swept exhaustively, cannot zero ANY '
      'failing check.' if best[0] == 0 else '')
print('The only bits with real leverage are the 3 structural controls, and those '
      'are exactly the branch flips already measured (x_2081, x_4287, x_24601: '
      '106 / 34 / 83 failing in the witness frame).')
