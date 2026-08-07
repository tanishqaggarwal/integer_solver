"""Solve the lattice over the CORRECT 12 coordinates, then test minimal representatives.

Coordinates (measured live in dimcheck.py, not assumed):
   live cofactors  1329, 9413, 10903, 17325     (free in the base harness)
   broken-atom wires 642, 28730, 29854, 31864   (NOT free in the base harness; assignable
                                                 here only because the deliverable breaks
                                                 their defining atoms -- engine3's demotion)
   carriers        7068, 4432, 9118, 8731
Far side is 12, not 13; the deliverable fails 7; the gap the cofactors buy is 5.

The magnitude question, stated so it can be answered: solve the target, then replace each
coordinate's delta by its balanced representative modulo the blocking modulus O gives, and
re-measure.  If collateral is magnitude-driven, the reduced solution scores better.
"""
import sys, os, json, time, collections, itertools
os.chdir('/home/user/integer_solver/solve_lab/agentM_work')
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentE_work')
sys.set_int_max_str_digits(20_000_000)
import harness as H
import engine as EB, engine3 as E3
import price as PR, fscore, sparse, sweep

DC = json.load(open('dimcheck.json'))
COORDS = DC['coords']
FAR12 = DC['far_side_failures']
D0 = json.load(open('DELTA0_FOR_M.json'))
p = int(D0['p'])
M642 = 2458959
M7068 = 7376877                      # 3 * 2458959

MOD = {642: M642, 17325: M642 * p, 1329: p, 9413: p, 10903: p, 7068: M7068}

vd = PR.load_deliverable()
freed, demote = PR.closure([642, 28730, 29854, 31864])
eng = E3.Eng(demote)
pinset = set(eng.pin)
seed = {f: vd[f] for f in eng.FREE if vd[f] != 0}
v0 = eng.forward(seed)
bad0 = eng.badatoms(v0)
F7 = sorted(fscore.fails(bad0))
print(f'deliverable: {fscore.score(bad0)} / failing {F7}', flush=True)
print(f'coordinates ({len(COORDS)}): {COORDS}', flush=True)

t0 = time.time()
cols = {}; aff = []
for f in COORDS:
    o = v0[f]
    b1, _ = sweep.fast_resid(eng, v0, bad0, {f: o + 1}, pinset)
    b2, _ = sweep.fast_resid(eng, v0, bad0, {f: o + 2}, pinset)
    b7, _ = sweep.fast_resid(eng, v0, bad0, {f: o + 7}, pinset)
    col = {}; ok = True
    for a in set(b1) | set(b2) | set(b7) | set(bad0):
        d1 = b1.get(a, 0) - bad0.get(a, 0)
        if b2.get(a, 0) - bad0.get(a, 0) != 2 * d1 or b7.get(a, 0) - bad0.get(a, 0) != 7 * d1:
            ok = False; break
        if d1:
            col[a] = d1
    if ok:
        aff.append(f); cols[f] = col
print(f'affine coordinates: {len(aff)}/{len(COORDS)} ({time.time()-t0:.0f}s)', flush=True)
for f in aff:
    print(f'   x_{f:<6d} -> atoms {sorted(cols[f])}', flush=True)
nonaff = [f for f in COORDS if f not in cols]
if nonaff:
    print(f'   NOT affine (excluded from the linear solve): {nonaff}', flush=True)


def rowfor(e):
    cm = collections.defaultdict(int); const = 0
    for c, a in H.eqt[e][2]:
        if a < 0:
            const += c
        else:
            cm[a] += c
    row = {}
    for f in aff:
        co = 0
        for a, d in cols[f].items():
            c = cm.get(a)
            if c:
                co += c * d
        if co:
            row[f] = co
    return row, -(const + sum(c * bad0[a] for a, c in cm.items() if a in bad0))


def measure(ch, label, dump=False):
    s = dict(seed)
    for k, v in ch.items():
        if v:
            s[k] = v
        else:
            s.pop(k, None)
    v = eng.forward(s)
    av = eng.badatoms(v)
    fl = sorted(fscore.fails(av))
    sc = PR.NEQ - len(fl)
    bits = {k: (v0[k] - ch[k]).bit_length() for k in ch}
    print(f'  [{label}] score {sc}  failing {len(fl)}  nbad {len(av)}  '
          f'max shift {max(bits.values()) if bits else 0} bits', flush=True)
    if sc > 39026:
        fn = f'M_min_{label}_{sc}.json'
        json.dump({f"x_{i}": int(v[i]) for i in range(PR.NV) if v[i] != 0}, open(fn, 'w'))
        print(f'      *** ABOVE BASELINE -> {fn} ***', flush=True)
    return sc, fl


def bal(x, m):
    r = x % m
    return r - m if r > m // 2 else r


print('\n=== solve the 7 failing equations over the 12 coordinates ===', flush=True)
rows = []; rhs = []
for e in F7:
    r, b = rowfor(e)
    rows.append(r); rhs.append(b)
idx = [i for i in range(len(rows)) if rows[i] or rhs[i] != 0]
sol, msg, _ = sparse.solve_sparse([rows[i] for i in idx], [rhs[i] for i in idx],
                                  verbose=False, maxcore=600, maxcorebits=8_000_000)
if sol is None:
    print(f'  NO SOLUTION ({msg})', flush=True)
else:
    ch = {f: v0[f] + d for f, d in sol.items() if d}
    print(f'  solved; {len(ch)} coordinates move, deltas (bits): '
          f'{ {f: (ch[f]-v0[f]).bit_length() for f in ch} }', flush=True)
    measure(ch, 'raw')

    print('\n=== the same solution with MINIMAL REPRESENTATIVES ===', flush=True)
    chm = dict(ch)
    for f in list(chm):
        m = MOD.get(f, p)
        d = chm[f] - v0[f]
        dr = bal(d, m)
        chm[f] = v0[f] + dr
    print(f'  reduced deltas (bits): { {f: (chm[f]-v0[f]).bit_length() for f in chm} }',
          flush=True)
    measure(chm, 'reduced')

    print('\n=== reduce ONE coordinate at a time (isolate which magnitude matters) ===',
          flush=True)
    for f in sorted(ch):
        c2 = dict(ch)
        m = MOD.get(f, p)
        c2[f] = v0[f] + bal(c2[f] - v0[f], m)
        measure(c2, f'only_x{f}')

print('\nbaseline 39026 ; far side 39021 (12 failing)', flush=True)
