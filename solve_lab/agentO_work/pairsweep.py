"""Agent O: sweep (a-bit, b-bit) selector PAIRS drawn from the 106 proven pin-solvable
   bits through E's full11 (1,1) simultaneous solve.  Writes jsonl to agentO_work/runs.

   usage: pairsweep.py <out.jsonl> <pairs...>   with pair = "a,b"
"""
import sys, os, json, time
ED = '/home/user/integer_solver/solve_lab/agentE_work'
OD = '/home/user/integer_solver/solve_lab/agentO_work'
sys.path.insert(0, ED)
os.chdir(ED)
sys.set_int_max_str_digits(20_000_000)
import engine as E, full11 as F

out = sys.argv[1]
if not os.path.isabs(out):
    out = OD + '/' + out
pairs = [tuple(int(x) for x in p.split(',')) for p in sys.argv[2:]]
BESTFILE = os.environ.get('OBEST', OD + '/best_pair')
THRESH = int(os.environ.get('OTHRESH', '7'))

f = open(out, 'a', buffering=1)
for a, b in pairs:
    t0 = time.time()
    try:
        r = F.solve_pair(a, b, verbose=False)
    except Exception as e:
        rec = {'a': a, 'b': b, 'err': f'{type(e).__name__}: {e}'[:120], 't': round(time.time() - t0, 1)}
        f.write(json.dumps(rec) + '\n'); print(rec, flush=True); continue
    if r is None:
        rec = {'a': a, 'b': b, 'err': 'nosol', 't': round(time.time() - t0, 1)}
        f.write(json.dumps(rec) + '\n'); print(rec, flush=True); continue
    n, ns, av = r
    rec = {'a': a, 'b': b, 'fails': n, 'score': 39033 - n, 'bad': av, 't': round(time.time() - t0, 1)}
    f.write(json.dumps(rec) + '\n'); print(rec, flush=True)
    if n < THRESH:
        v = E.forward(ns)
        tag = f'{BESTFILE}_{a}_{b}_{39033-n}'
        json.dump({f"x_{i}": str(int(v[i])) for i in range(E.NV) if v[i] != 0}, open(tag + '.json', 'w'))
        json.dump({str(k): str(int(x)) for k, x in ns.items()}, open(tag + '_seed.json', 'w'))
        print('*** WROTE', tag, flush=True)
