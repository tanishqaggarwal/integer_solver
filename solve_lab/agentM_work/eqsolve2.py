"""Equation-space solve with an ITERATIVELY WIDENED knob set.

eqsolve.py showed the objective is right but the 19-knob set leaves 7 equations
unsatisfiable, with a divisibility obstruction on eq 29125.  Widen the knob set by
expanding cones round by round and re-solve, to see whether the obstruction is an
artifact of the knob set or intrinsic.
"""
import sys, os, json, time, pickle, collections
os.chdir('/home/user/integer_solver/solve_lab/agentM_work')
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentE_work')
sys.set_int_max_str_digits(20_000_000)
import harness as H
import engine as E_broken
import engine2 as E2, fast2, mcore2 as M, chan2 as C, sparse

NEQ = len(H.eqt)
ROUNDS = int(sys.argv[1]) if len(sys.argv) > 1 else 2
MAXK = int(sys.argv[2]) if len(sys.argv) > 2 else 1200

vd = M.load_vec()
base = E2.seed_of(vd)
v0 = E2.forward(base)
bad0 = E2.badatoms(v0)
fails0 = set(E2.eqfails(bad0))
print(f'base score {NEQ-len(fails0)}, failing {sorted(fails0)}', flush=True)

atom2eq = collections.defaultdict(list)
for e, (issq, outer, terms) in enumerate(H.eqt):
    for c, a in terms:
        if a >= 0:
            atom2eq[a].append(e)

FREESET = set(E2.FREE)
cand = set(E2.PIN) | set(C.CLUSTERKN)
for a in sorted(bad0):
    cand |= set(E_broken.cone(a)[1])
# widen: cones of every atom appearing in the currently-failing equations
for e in sorted(fails0):
    issq, outer, terms = H.eqt[e]
    for c, a in terms:
        if a >= 0:
            cand |= set(E_broken.cone(a)[1])
cand = sorted(f for f in cand if f in FREESET)
print(f'round-0 candidates: {len(cand)}', flush=True)

allaff = {}
seen = set()
for rnd in range(ROUNDS + 1):
    newc = [f for f in cand if f not in seen]
    if not newc:
        break
    seen |= set(newc)
    t0 = time.time()
    aff, cols = C.affine_cols(v0, bad0, newc)
    allaff.update(cols)
    print(f'  round {rnd}: probed {len(newc)}, affine {len(aff)}, total knobs {len(allaff)} ({time.time()-t0:.0f}s)', flush=True)
    if len(allaff) > MAXK:
        break
    # expand: atoms these knobs touch -> their cones
    ext = set()
    for f in aff:
        for a in cols[f]:
            ext |= set(E_broken.cone(a)[1])
    cand = sorted(f for f in ext if f in FREESET and f not in seen)
    if len(seen) + len(cand) > 4000:
        cand = cand[:4000 - len(seen)]

knobs = sorted(allaff)
print(f'\n{len(knobs)} affine knobs total', flush=True)

touched = set(bad0)
for f in knobs:
    touched |= set(allaff[f])
eqs = set(fails0)
for a in touched:
    eqs |= set(atom2eq[a])
eqs = sorted(eqs)
print(f'{len(touched)} touched atoms, {len(eqs)} equations', flush=True)

def eqval(av, e):
    issq, outer, terms = H.eqt[e]
    s = 0
    for c, a in terms:
        if a < 0:
            s += c
        elif a in av:
            s += c * av[a]
    return s

rows, rhs = [], []
for e in eqs:
    issq, outer, terms = H.eqt[e]
    cmap = {}
    for c, a in terms:
        if a >= 0:
            cmap[a] = cmap.get(a, 0) + c
    row = {}
    for f in knobs:
        coef = 0
        for a, d in allaff[f].items():
            c = cmap.get(a)
            if c:
                coef += c * d
        if coef:
            row[f] = coef
    rows.append(row); rhs.append(-eqval(bad0, e))

nodep = [eqs[i] for i in range(len(eqs)) if not rows[i] and rhs[i] != 0]
print(f'failing eqs with NO knob dependence: {len(nodep)} {nodep[:10]}', flush=True)

idx = [i for i in range(len(rows)) if rows[i]]
print(f'solving {len(idx)} rows...', flush=True)
sol, msg, _ = sparse.solve_sparse([rows[i] for i in idx], [rhs[i] for i in idx],
                                  names=[eqs[i] for i in idx], verbose=False,
                                  maxcore=600, maxcorebits=8_000_000)
print('full solve:', 'FOUND' if sol else f'none ({msg})', flush=True)

def tryapply(sol, tag):
    if not sol:
        return None
    ns = dict(base)
    for f, d in sol.items():
        if d:
            ns[f] = v0[f] + d
    v = E2.forward(ns)
    av = E2.badatoms(v)
    sc = NEQ - len(E2.eqfails(av))
    print(f'  [{tag}] score {sc} nbad {len(av)}', flush=True)
    if sc > 39026:
        C.dump(v, f'M_eq2_{sc}.json')
        print(f'  *** ABOVE BASELINE -> M_eq2_{sc}.json ***', flush=True)
    return sc

tryapply(sol, 'full')

# greedy over rows, prioritising the currently-failing equations first
order = sorted(idx, key=lambda i: (eqs[i] not in fails0))
keep = []
t0 = time.time()
for i in order:
    cand_idx = keep + [i]
    s2, _, _ = sparse.solve_sparse([rows[j] for j in cand_idx], [rhs[j] for j in cand_idx],
                                   verbose=False, maxcore=600, maxcorebits=8_000_000)
    if s2 is not None:
        keep = cand_idx
    if time.time() - t0 > 1500:
        print('  [greedy budget]', flush=True); break
kept_fail = [eqs[j] for j in keep if eqs[j] in fails0]
print(f'greedy kept {len(keep)}/{len(order)}; of the 7 failing it satisfies {len(kept_fail)}: {kept_fail}', flush=True)
if keep:
    s3, _, _ = sparse.solve_sparse([rows[j] for j in keep], [rhs[j] for j in keep],
                                   verbose=False, maxcore=600, maxcorebits=8_000_000)
    tryapply(s3, 'greedy')
print('done', flush=True)
