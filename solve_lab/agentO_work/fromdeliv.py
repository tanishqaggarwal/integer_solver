"""Start from the DELIVERABLE state (39,026) and run the full-knob lazy simultaneous
   repair on its 8 bad atoms.  The deliverable is reproduced from its free variables."""
import sys, os, json, time
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentO_work')
import simO, engine as E, harness as H, sparse, lazy
OD = '/home/user/integer_solver/solve_lab/agentO_work'

d = json.load(open('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json'))
vd = [0] * E.NV
for k, x in d.items():
    vd[int(k.split('_')[1])] = int(x)
FR = [u for u in range(E.NV) if E.definer[u] is None]
seed = {u: vd[u] for u in FR if vd[u] != 0}
v = E.forward(seed)
bad = E.badatoms(v)
print('reconstructed from', len(seed), 'free vars: bad =', sorted(bad),
      'fails =', len(E.eqfails(bad)), flush=True)
assert all(v[u] == vd[u] for u in range(E.NV)), 'forward does not reproduce deliverable'
print('forward reproduces the deliverable exactly', flush=True)

maxr = int(sys.argv[1]) if len(sys.argv) > 1 else 6
maxv = int(sys.argv[2]) if len(sys.argv) > 2 else 8000
frozen = set(int(x) for x in sys.argv[3].split(',')) if len(sys.argv) > 3 and sys.argv[3] else set()
t0 = time.time()
r = lazy.run(seed, frozen, maxr=maxr, maxv=maxv, iters=40)
print(f'BEST fails={r[0]} score={39033-r[0]} bad={r[2]} ({time.time()-t0:.0f}s)', flush=True)
if r[0] < 7:
    json.dump({f"x_{i}": str(int(r[3][i])) for i in range(E.NV) if r[3][i] != 0},
              open(f'{OD}/deliv_repair_{39033-r[0]}.json', 'w'))
    print('*** WROTE improvement', flush=True)
