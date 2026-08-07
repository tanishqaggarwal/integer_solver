"""THE GATE: does forward from the deliverable's free inputs reproduce 39,026
with its 8 nonzero atoms intact?"""
import sys, os, json
os.chdir('/home/user/integer_solver/solve_lab/agentM_work')
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentE_work')
sys.set_int_max_str_digits(20_000_000)
import harness as H
import engine2 as E2

BAD8 = {23616, 23617, 36659, 36660, 36661, 36662, 36663, 36664}
DELIV = '/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json'


def load_vec(path):
    d = json.load(open(path))
    v = [0] * H.NV
    for k, val in d.items():
        v[int(k.split('_')[1])] = int(val)
    return v


vd = load_vec(DELIV)
print('PIN (vars promoted to free):', E2.PIN)
print('|FREE| E =', len(set(H.FREE)), ' -> engine2 FREE =', len(E2.FREE),
      ' |SEQ| ', len(H.SEQ), '->', len(E2.SEQ))

seed = E2.seed_of(vd)
print('seed size (nonzero free inputs):', len(seed))

v = E2.forward(seed)
diff = [u for u in range(H.NV) if v[u] != vd[u]]
print()
print('=== GATE 1: exact vector reproduction ===')
print('vars differing from deliverable:', len(diff), diff[:20])

av = E2.badatoms(v)
fails = E2.eqfails(av)
print()
print('=== GATE 2: score ===')
print(f'satisfied {len(H.eqt) - len(fails)}/{len(H.eqt)}   (failing {len(fails)})')
print('failing eqs:', sorted(fails))

print()
print('=== GATE 3: the 8 nonzero atoms intact ===')
nz = set(av)
print('nonzero atom count:', len(nz))
print('nonzero atoms      :', sorted(nz))
print('expected BAD8      :', sorted(BAD8))
print('MATCH:', nz == BAD8)

ok = (len(diff) == 0 and len(H.eqt) - len(fails) == 39026 and nz == BAD8)
print()
print('GATE PASSED' if ok else 'GATE FAILED')
