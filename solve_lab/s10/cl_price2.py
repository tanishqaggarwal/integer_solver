"""CL: minimum equation sacrifice with LEGAL columns only (boolean-pinned inputs frozen).
Atom-level closure -> upper bound on the reachable score from this frame."""
import os, sys, json, time, collections
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad
import cl_engine as E
P = E.P
definer, atom_out = L.definer, L.atom_out
v0 = L.load(os.path.join(HERE,'mod9118_0.json'))
av0 = L.all_atom_values(v0)
BADA = [a for a in range(L.NA) if a not in atom_out and av0[a]]
GAIN = set()
for a in BADA: GAIN |= set(L.atom2eq.get(a, {}))

BOOLV = set()
for a in range(L.NA):
    if a in atom_out: continue
    vs = set()
    for mo in L.polys[a]: vs |= set(mo)
    if len(vs) == 1 and max(len(mo) for mo in L.polys[a]) == 2:
        BOOLV.add(next(iter(vs)))

D = json.load(open(os.path.join(HERE, 'cl_closure_cols.json')))
cols = {int(u): {int(c): int(d) for c, d in mm.items()} for u, mm in D['cols'].items()}
allrows = sorted(set().union(*[set(c) for c in cols.values()]) | set(BADA))

def run(U, tag, steps=26):
    m = len(U)
    print(f'\n===== {tag}: {len(allrows)} atom rows x {m} columns =====')
    def analyse(drop, track=True):
        rows = [a for a in allrows if a not in drop]
        n = len(rows); ri = {a: i for i, a in enumerate(rows)}
        A = [[0]*(m+1) for _ in range(n)]
        for j, u in enumerate(U):
            for a, d in cols[u].items():
                if a in ri: A[ri[a]][j] = d % P
        for a in rows: A[ri[a]][m] = (-av0[a]) % P
        Tr = [{i: 1} for i in range(n)] if track else None
        r = 0
        for c in range(m):
            pr = next((i for i in range(r, n) if A[i][c]), None)
            if pr is None: continue
            A[r], A[pr] = A[pr], A[r]
            if Tr: Tr[r], Tr[pr] = Tr[pr], Tr[r]
            inv = pow(A[r][c], -1, P)
            if inv != 1:
                A[r] = [x*inv % P for x in A[r]]
                if Tr: Tr[r] = {k: x*inv % P for k, x in Tr[r].items()}
            Ar = A[r]
            for i in range(r+1, n):
                f = A[i][c]
                if not f: continue
                A[i] = [(A[i][k]-f*Ar[k]) % P for k in range(m+1)]
                if Tr:
                    ti = Tr[i]
                    for k, x in Tr[r].items():
                        nv = (ti.get(k, 0)-f*x) % P
                        if nv: ti[k] = nv
                        elif k in ti: del ti[k]
            r += 1
        obs = [i for i in range(r, n) if A[i][m] % P and not any(A[i][j] for j in range(m))]
        certs = [{rows[k] for k, w in Tr[i].items() if w % P} for i in obs] if Tr else []
        return r, obs, certs
    drop = set()
    for step in range(steps):
        r, obs, certs = analyse(drop)
        lost = set()
        for a in drop: lost |= set(L.atom2eq.get(a, {}))
        fail = len(lost | (GAIN if obs else set()))
        print(f'  step {step}: rank {r}, obstructions {len(obs)}, sacrificed atoms {len(drop)}, '
              f'equations lost {len(lost)}  -> if solved, score {L.NEQ - len(lost)}')
        if not obs:
            print(f'  *** CONSISTENT.  sacrifice {sorted(drop)} costs {len(lost)} equations; '
                  f'best reachable score {L.NEQ - len(lost)}')
            json.dump({'drop': sorted(drop), 'eqs_lost': sorted(lost)},
                      open(os.path.join(HERE, f'cl_sac_{m}.json'), 'w'))
            return L.NEQ - len(lost)
        pool = collections.Counter()
        for cc in certs: pool.update(cc)
        opts = [a for a in pool if a not in BADA and a not in drop]
        if not opts:
            print('  no further sacrifice candidates'); return None
        opts.sort(key=lambda a: (len(set(L.atom2eq.get(a, {})) - lost - GAIN), -pool[a]))
        best = None
        for a in opts[:8]:
            r2, obs2, _ = analyse(drop | {a}, track=False)
            l2 = set(lost) | set(L.atom2eq.get(a, {}))
            if best is None or (len(obs2), len(l2)) < best[0]:
                best = ((len(obs2), len(l2)), a)
        drop.add(best[1])
    return None

UALL = sorted(cols)
ULEG = [u for u in UALL if u not in BOOLV]
run(ULEG, f'LEGAL columns ({len(ULEG)})')
run(UALL, f'ALL columns ({len(UALL)})')
