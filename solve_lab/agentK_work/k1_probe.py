#!/usr/bin/env python3
"""K1: forward-evaluate the deliverable, read off the root stage, test the chordK law."""
import sys, os, json, pickle, time
F = '/home/user/integer_solver/solve_lab/agentF_work'
sys.path.insert(0, F)
from fwd import Engine, NV

p = 115792089237316195423570985008687907853269984665640564039457584007908834671663
K = 97553848499418123410591666447050222001188385549510401465815187079080512838891

t0 = time.time()
E = Engine()
print('engine built', time.time() - t0, 'free', len(E.free), 'res', len(E.res))

d = json.load(open('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json'))
full = [0] * NV
for k, v in d.items():
    full[int(k[2:])] = int(v)

# take only free inputs from the deliverable, forward-evaluate the rest
v = [0] * NV
freeset = set(E.free)
for i in E.free:
    v[i] = full[i]
r = E.run(v)
bad = E.score(r)
print('fwd from deliverable free inputs: failing eqs', len(bad), '=> score', 39033 - len(bad))
nz = [i for i, x in enumerate(r) if x]
print('nonzero residual atoms:', len(nz), [E.res[i] for i in nz][:20])

# agreement with deliverable on defined vars
diff = [i for i in range(NV) if v[i] != full[i]]
print('vars differing from deliverable:', len(diff), diff[:20])

json.dump({'x_%d' % i: v[i] for i in range(NV) if v[i]}, open('/home/user/integer_solver/solve_lab/agentK_work/fwd_deliverable.json', 'w'))
pickle.dump(v, open('/home/user/integer_solver/solve_lab/agentK_work/v_deliverable.pkl', 'wb'))

roles = json.load(open(F + '/stage_roles.json'))
def chordK(ax, ay, bx, by):
    l = (by - ay) * pow(bx - ax, p - 2, p) % p
    ox = (l * l - ax - bx - K) % p
    oy = (l * (ax - ox) - ay) % p
    return ox, oy

r0 = roles['15298'][0]
ax, ay = v[r0['inA'][0]], v[r0['inA'][1]]
bx, by = v[r0['inB'][0]], v[r0['inB'][1]]
ox, oy = v[r0['out'][0]], v[r0['out'][1]]
print('root inA', ax % p, ay % p)
print('root inB', bx % p, by % p)
print('root out', ox % p, oy % p)
for (A, B) in [((ax, ay), (bx, by)), ((ay, ax), (by, bx))]:
    try:
        c = chordK(A[0] % p, A[1] % p, B[0] % p, B[1] % p)
        print('chordK ->', c, 'matches out(x,y)?', c == (ox % p, oy % p), 'matches out(y,x)?', c == (oy % p, ox % p))
    except Exception as e:
        print('err', e)
print('done', time.time() - t0)
