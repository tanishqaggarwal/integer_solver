"""S10 step 99: beam search over the ENRICHED move set.

Moves:
  (Z) zero-atom, realised directly or through a gate's free input (2-level)
  (N) mod-p Newton on a free input  (fix a residue, not a value)
  (H) p-handle absorption after (N)
Score = equations satisfied.  Negative steps allowed; that is the whole point --
every previous search was a pure hill-climb and every wall it hit was a plateau.
"""
import os, sys, time, heapq, json
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
P = ad.P
definer, atom_out = L.definer, L.atom_out
FREE = set(ad.FREE)
FORBID = {2081, 4287}
DEADLINE = time.time() + float(os.environ.get('BUDGET', '3000'))

def score(v):
    return L.NEQ - len(L.failing_eqs(L.all_atom_values(v)))

def zero_moves(a, v):
    out = []
    for w in sorted(set(L.avars[a])):
        if w in FORBID: continue
        tgt = T.solve_lin(a, w, v)
        if tgt is None or tgt == v[w]: continue
        if w in FREE:
            out.append((w, tgt))
        else:
            d = definer.get(w)
            if d is None: continue
            vv = list(v); vv[w] = tgt
            for u in sorted(set(L.avars[d])):
                if u == w or u not in FREE or u in FORBID: continue
                nv = T.solve_lin(d, u, vv)
                if nv is not None: out.append((u, nv))
    return out

def newton_moves(a, v, av, topn=10):
    vm = [x % P for x in v]
    try: g = ad.grad(a, vm)
    except Exception: return []
    r = av[a] % P
    if r == 0: return []
    cand = []
    for u, d in g.items():
        if u in FORBID or d % P == 0: continue
        cand.append((len(L.var_atoms[u]), u, (-r * pow(d, -1, P)) % P))
    cand.sort()
    return [(u, v[u] + delta) for _, u, delta in cand[:topn]]

def expand(v):
    av = L.all_atom_values(v)
    nz = [a for a in range(L.NA) if av[a]]
    seen, kids = set(), []
    for a in nz:
        for u, nv in zero_moves(a, v) + newton_moves(a, v, av):
            if (u, nv) in seen: continue
            seen.add((u, nv))
            w = list(v); w[u] = nv
            ad.fwd(w, rounds=6)
            # opportunistic handle absorption on every still-nonzero atom
            for _ in range(3):
                aw = L.all_atom_values(w)
                did = False
                for b in [c for c in range(L.NA) if aw[c]]:
                    for u2, nv2 in zero_moves(b, w):
                        t2 = list(w); t2[u2] = nv2
                        ad.fwd(t2, rounds=6)
                        if score(t2) > score(w):
                            w = t2; did = True; break
                    if did: break
                if not did: break
            kids.append((score(w), (u, nv), w))
    return kids

if __name__ == '__main__':
    start = L.load(os.path.join(HERE, sys.argv[1] if len(sys.argv) > 1 else 'mod9118_0.json'))
    s0 = score(start)
    print(f'seed score {s0}', flush=True)
    BEAM = 5
    frontier = [(s0, start, [])]
    best = (s0, start, [])
    visited = set()
    depth = 0
    while frontier and time.time() < DEADLINE and depth < 8:
        depth += 1
        kids = []
        for s, v, path in frontier:
            if time.time() > DEADLINE: break
            for ks, mv, w in expand(v):
                key = (ks, tuple(sorted((mv[0],))), len(path))
                kids.append((ks, w, path + [mv]))
                if ks > best[0]:
                    best = (ks, w, path + [mv])
                    T.save(w, os.path.join(HERE, f'beam_{ks}.json'))
                    print(f'  depth {depth}: NEW BEST {ks}  path {[m[0] for m in best[2]]}',
                          flush=True)
        if not kids: break
        kids.sort(key=lambda t: -t[0])
        # de-duplicate by resulting nonzero-atom signature
        frontier, sigs = [], set()
        for ks, w, path in kids:
            aw = L.all_atom_values(w)
            sig = tuple(b for b in range(L.NA) if aw[b])
            if sig in sigs: continue
            sigs.add(sig); frontier.append((ks, w, path))
            if len(frontier) >= BEAM: break
        print(f'depth {depth}: frontier {[f[0] for f in frontier]}  best {best[0]}',
              flush=True)
    print(f'FINAL best {best[0]} via {[m[0] for m in best[2]]}')
