"""S10 step 26: wide beam search -- is ANY d != 0 (mod p) closable for x_28730?

If yes, congruence (2) dies and the score is 39,027 (and the endgame changes:
only one congruence would remain).  d = p closes trivially but does not move K2.
"""
import os, sys, json, time, heapq
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
FORBID = {28730, 4432}   # do not let the search simply undo the seed move


def bad_after(v, changed):
    touched = L.touched_atoms(v, base_av, changed)
    return sorted(set(a for a, val in touched.items() if val != 0) - BASE_NZ)


def search(d, tag):
    v = list(base)
    ch, _ = L.ripple(v, {28730: base[28730] + d, 4432: base[4432] + d}, block=BLOCK)
    bad = bad_after(v, ch)
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
                    key = tuple(b2)
                    if not b2:
                        av = L.all_atom_values(w)
                        fail = L.failing_eqs(av)
                        print(f'  *** {tag} CLOSED at depth {depth} via x_{u}: '
                              f'nz={[x for x in range(L.NA) if av[x]]} '
                              f'failing={len(fail)} score={L.NEQ-len(fail)} '
                              f'K2={w[28730] % P}', flush=True)
                        T.save(w, os.path.join(HERE, f'freed28730_{tag}.json'))
                        return w
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
for d, tag in ((1, 'd1'), (-1, 'dm1'), (2, 'd2'), (7376877, 'd7376877'),
               (P + 1, 'dp1'), (12846437, 'd12846437')):
    print(f'\n===== d = {d if abs(d) < 10**9 else "p+1"} =====', flush=True)
    r = search(d, tag)
    print(f'  -> {"CLOSED" if r is not None else "not closed"}  ({time.time()-t0:.0f}s)',
          flush=True)
