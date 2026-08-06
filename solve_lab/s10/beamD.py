"""S10 step 27: same wide beam search for congruence (1) -- move x_7068 by d not= 0 (mod p).

d = k*p is already known free (s10/repairD.py).  The question is whether any
d with d % p != 0 can be absorbed, which would kill congruence (1).
"""
import os, sys, json, time
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L
import tools as T

P = 2**256 - 2**32 - 977
NZ = [22229, 22230, 35758, 35759, 35760, 35761, 35762]
BLOCK = set(NZ) | {22231}
base = L.load(os.path.join(LAB, 'best', 'new_instance_partial_39026.json'))
base_av = L.all_atom_values(base)
BASE_NZ = set(a for a in range(L.NA) if base_av[a])
BEAM = int(os.environ.get('BEAM', 120))
DEPTH = int(os.environ.get('DEPTH', 9))


def bad_after(v, changed):
    touched = L.touched_atoms(v, base_av, changed)
    return sorted(set(a for a, val in touched.items() if val != 0) - BASE_NZ)


def search(seeds, tag, forbid):
    global FORBID
    FORBID = set(forbid)
    v = list(base)
    ch, _ = L.ripple(v, dict(seeds), block=BLOCK)
    bad = bad_after(v, ch)
    print(f'  {tag} initial collateral {bad}', flush=True)
    frontier = [(v, ch, bad, ())]
    seen = set()
    for depth in range(DEPTH):
        nxt = []
        for vv, cc, bb, path in frontier:
            for a in bb:
                for u in sorted(L.avars[a]):
                    if (a, u) in path or u in FORBID:
                        continue
                    nv = T.solve_lin(a, u, vv)
                    if nv is None or nv == vv[u]:
                        continue
                    w = list(vv)
                    try:
                        ch2, _ = L.ripple(w, {u: nv}, block=BLOCK)
                    except Exception:
                        continue
                    allch = dict(cc); allch.update(ch2)
                    b2 = bad_after(w, allch)
                    if not b2:
                        av = L.all_atom_values(w)
                        fail = L.failing_eqs(av)
                        D = w[7068] - w[2099]
                        print(f'  *** {tag} CLOSED depth {depth} via x_{u}: '
                              f'nz={[x for x in range(L.NA) if av[x]]} '
                              f'failing={len(fail)} score={L.NEQ-len(fail)} '
                              f'D%p={D % P}', flush=True)
                        T.save(w, os.path.join(HERE, f'freedD_{tag}.json'))
                        return w
                    key = tuple(b2)
                    if key in seen:
                        continue
                    seen.add(key)
                    nxt.append((w, allch, b2, path + ((a, u),)))
        nxt.sort(key=lambda t: (len(t[2]), sum(len(L.atom2eq.get(a, {})) for a in t[2])))
        frontier = nxt[:BEAM]
        if not frontier:
            print(f'  {tag} depth {depth}: exhausted', flush=True)
            return None
        print(f'  {tag} depth {depth}: {len(nxt)} branches, best={frontier[0][2]}',
              flush=True)
    return None


t0 = time.time()
FORBID = set()
CASES = [
    ({7068: base[7068] + 1}, 'x7068+1'),
    ({7068: base[7068] - 1}, 'x7068-1'),
    ({7068: base[7068] + 12846437}, 'x7068+12846437'),
    # move the mirror together with it -- a29539 wants x_1308 == x_14853 (mod p)
    ({7068: base[7068] + 1, 14853: base[14853] + 1}, 'x7068+1,x14853+1'),
    ({2099: base[2099] + 1}, 'x2099+1'),
    ({37158: base[37158] + 1}, 'x37158+1'),
    ({10878: base[10878] + 1}, 'x10878+1'),
    ({22542: base[22542] + 1}, 'x22542+1'),
]
for seeds, tag in CASES:
    print(f'\n===== {tag} =====', flush=True)
    r = search(seeds, tag, set(seeds) | {7068, 2099})
    print(f'  -> {"CLOSED" if r is not None else "not closed"} ({time.time()-t0:.0f}s)',
          flush=True)
