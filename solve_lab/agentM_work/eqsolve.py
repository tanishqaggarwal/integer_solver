"""Solve in EQUATION space, not atom space.

simsolve targets av[a]=0 for every bad atom -- that is why it loses (39,008): the
39,026 deliverable keeps 8 atoms NONZERO precisely because they cancel inside the
equations.  Correct objective: choose knob deltas so that each equation's total is
zero, allowing nonzero atoms.  The 5 vars engine2 promotes to free are the natural
knobs -- they drive 7 of the 8 residual atoms affinely.
"""
import sys, os, json, time, pickle, collections
os.chdir('/home/user/integer_solver/solve_lab/agentM_work')
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentE_work')
sys.set_int_max_str_digits(20_000_000)
import harness as H
import engine2 as E2, fast2, mcore2 as M, chan2 as C, sparse

NEQ = len(H.eqt)
vd = M.load_vec()
base = E2.seed_of(vd)
v0 = E2.forward(base)
bad0 = E2.badatoms(v0)
fails0 = E2.eqfails(bad0)
print(f'base: score {NEQ-len(fails0)}  bad {sorted(bad0)}  failing {sorted(fails0)}', flush=True)

# atom -> equations index
atom2eq = collections.defaultdict(list)
for e, (issq, outer, terms) in enumerate(H.eqt):
    for c, a in terms:
        if a >= 0:
            atom2eq[a].append(e)

print('\nequations touched by each bad atom:')
for a in sorted(bad0):
    es = atom2eq[a]
    bad_es = [e for e in es if e in set(fails0)]
    print(f'  atom {a}: {len(es)} eqs, of which failing: {bad_es}')

# ---- candidate knobs: the 5 promoted vars + cluster knobs + cone of the bad atoms ----
import engine as E_broken
cand = set(E2.PIN) | set(C.CLUSTERKN)
for a in sorted(bad0):
    cand |= set(E_broken.cone(a)[1])
cand = sorted(f for f in cand if (f in set(E2.FREE)) and f not in ())
print(f'\n{len(cand)} candidate knobs', flush=True)

t0 = time.time()
aff, cols = C.affine_cols(v0, bad0, cand)
print(f'affine knobs: {len(aff)}  ({time.time()-t0:.0f}s)', flush=True)
print('  PIN knobs affine?', {f: (f in cols) for f in E2.PIN})

# ---- build equation-space system ----
# equations that can change: those touching any atom in any knob column, plus current failures
touched_atoms = set(bad0)
for f in aff:
    touched_atoms |= set(cols[f])
eqs = set(fails0)
for a in touched_atoms:
    eqs |= set(atom2eq[a])
eqs = sorted(eqs)
print(f'{len(touched_atoms)} touched atoms, {len(eqs)} equations in the system', flush=True)

# current value of each equation
def eqval(av, e):
    issq, outer, terms = H.eqt[e]
    s = 0
    for c, a in terms:
        if a < 0:
            s += c
        elif a in av:
            s += c * av[a]
    return s

rhs = []
rows = []
for e in eqs:
    s0 = eqval(bad0, e)
    row = {}
    issq, outer, terms = H.eqt[e]
    cmap = {a: c for c, a in terms if a >= 0}
    for f in aff:
        coef = 0
        for a, d in cols[f].items():
            c = cmap.get(a)
            if c:
                coef += c * d
        if coef:
            row[f] = coef
    rows.append(row)
    rhs.append(-s0)

nz = sum(1 for r in rows if r)
print(f'rows with any knob dependence: {nz}/{len(rows)}', flush=True)
# drop rows that are all-zero AND already satisfied (nothing to do)
keepidx = [i for i, r in enumerate(rows) if r or rhs[i] != 0]
rows2 = [rows[i] for i in keepidx]
rhs2 = [rhs[i] for i in keepidx]
eqs2 = [eqs[i] for i in keepidx]
unfixable = [eqs[i] for i in keepidx if not rows[i] and rhs[i] != 0]
print(f'{len(rows2)} rows kept; {len(unfixable)} failing eqs have NO knob dependence: {unfixable[:20]}', flush=True)

solvable = [i for i in range(len(rows2)) if rows2[i]]
print(f'attempting exact solve on {len(solvable)} solvable rows', flush=True)
sol, msg, _ = sparse.solve_sparse([rows2[i] for i in solvable], [rhs2[i] for i in solvable],
                                  names=[eqs2[i] for i in solvable], verbose=False,
                                  maxcore=400, maxcorebits=5_000_000)
print('full solve:', 'FOUND' if sol else f'none ({msg})', flush=True)

best = (NEQ - len(fails0), None)
def tryapply(sol, tag):
    global best
    if not sol:
        return
    ns = dict(base)
    for f, d in sol.items():
        if d:
            ns[f] = v0[f] + d
    v = E2.forward(ns)
    av = E2.badatoms(v)
    sc = NEQ - len(E2.eqfails(av))
    print(f'  [{tag}] score {sc}  nbad {len(av)}', flush=True)
    if sc > best[0]:
        best = (sc, ns)
        C.dump(v, f'M_eq_{sc}.json')
        print(f'  *** NEW BEST {sc} -> M_eq_{sc}.json ***', flush=True)
    return sc

tryapply(sol, 'full')

# greedy maximal satisfiable subset
keep = []
t0 = time.time()
for i in solvable:
    idx = keep + [i]
    s2, _, _ = sparse.solve_sparse([rows2[j] for j in idx], [rhs2[j] for j in idx],
                                   verbose=False, maxcore=400, maxcorebits=5_000_000)
    if s2 is not None:
        keep = idx
    if time.time() - t0 > 1200:
        print('  [greedy time budget]', flush=True); break
print(f'greedy kept {len(keep)}/{len(solvable)} rows', flush=True)
if keep:
    s3, _, _ = sparse.solve_sparse([rows2[j] for j in keep], [rhs2[j] for j in keep],
                                   verbose=False, maxcore=400, maxcorebits=5_000_000)
    tryapply(s3, 'greedy')

print('\nBEST:', best[0], flush=True)
