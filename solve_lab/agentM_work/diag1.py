"""Step 1 diagnostic: why does E's forward from the deliverable's free inputs give 39008?"""
import sys, os, json, math, collections
os.chdir('/home/user/integer_solver/solve_lab/agentM_work')
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentE_work')
sys.set_int_max_str_digits(20_000_000)
import harness as H, engine as E

DELIV = '/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json'

def load_vec(path):
    d = json.load(open(path))
    v = [0] * H.NV
    for k, val in d.items():
        v[int(k.split('_')[1])] = int(val)
    return v

vd = load_vec(DELIV)
FREE = E.FREE
seed = {f: vd[f] for f in FREE if vd[f] != 0}
print('deliverable nonzero free inputs:', len(seed))

vf = E.forward(seed)
diff = [u for u in range(H.NV) if vf[u] != vd[u]]
print('vars differing after forward:', len(diff))
print('  all in SEQ (derived)?', all(H.definer[u] is not None for u in diff))
print('  diff vars:', diff)

# classify each differing var by its definer atom kind + root structure
posn = {u: k for k, u in enumerate(H.SEQ)}
diff_sorted = sorted(diff, key=lambda u: posn.get(u, -1))
print()
print('--- per-var analysis (in topo order) ---')
ns = {'v': list(vd), '__builtins__': {}}   # evaluate coefficients AT THE DELIVERABLE
vv = ns['v']
info = []
for u in diff_sorted:
    i, kind = H.definer[u]
    c = H.acodes[i]
    old = vv[u]
    vv[u] = 0; c0 = eval(c, ns)
    vv[u] = 1; c1 = eval(c, ns)
    vv[u] = 2; c2 = eval(c, ns)
    vv[u] = old
    A2 = c2 - 2 * c1 + c0
    if A2 == 0:
        sl = c1 - c0
        roots = [] if (sl == 0 or c0 % sl) else [-c0 // sl]
        k2 = 'lin'
    else:
        A = A2 // 2; B = c1 - c0 - A; C = c0
        disc = B * B - 4 * A * C
        k2 = 'quad'
        if disc < 0 or A == 0:
            roots = []
        else:
            r = math.isqrt(disc)
            roots = sorted({(-B + s) // (2 * A) for s in (r, -r)
                            if (-B + s) % (2 * A) == 0}) if r * r == disc else []
    info.append((u, i, kind[0], k2, len(roots), vd[u] in roots, vf[u], vd[u], roots))

nq = sum(1 for x in info if x[3] == 'quad')
n2 = sum(1 for x in info if x[4] == 2)
ndelivroot = sum(1 for x in info if x[5])
nfzero = sum(1 for x in info if x[6] == 0)
print(f'quad-defined: {nq}/{len(info)}   two-root: {n2}   deliverable value IS a root: {ndelivroot}   forward gave 0: {nfzero}')
print()
for u, i, kd, k2, nr, isroot, fv, dv, roots in info:
    rs = [str(r)[:24] + ('...' if len(str(r)) > 24 else '') for r in roots]
    print(f'x_{u:6d} atom{i:6d} kind={kd:4s}/{k2:4s} nroots={nr} delivIsRoot={isroot} '
          f'fwd={str(fv)[:14]:>14s} deliv={str(dv)[:14]:>14s} roots={rs}')

import pickle
pickle.dump({'diff': diff_sorted, 'info': info, 'seed': seed},
            open('diag1.pkl', 'wb'))
