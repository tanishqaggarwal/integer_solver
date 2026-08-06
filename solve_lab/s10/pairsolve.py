"""S11 step 5: solve BOTH cluster congruences exactly, on a pair of free inputs.

The residues a21617 and a29539 are exactly linear mod p in every non-boolean free
input (s10/linearity.py).  So for any pair (u1,u2) whose 2x2 gradient matrix is
invertible mod p there is an EXACT delta zeroing both congruences at once.  The
collateral checks are nonlinear, so the linear closure's veto does not apply --
settle them by construction instead.
"""
import os, sys, time, itertools, collections
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
P = ad.P
definer, atom_out = L.definer, L.atom_out
FREE = set(ad.FREE)
FORBID = {2081, 4287}
BAD = [21617, 29539]

def score(v):
    return L.NEQ - len(L.failing_eqs(L.all_atom_values(v)))

def settle(v, rounds=12):
    """apply every zero-collateral repair (solo handles) until nothing is left."""
    for _ in range(rounds):
        av = L.all_atom_values(v)
        nz = [a for a in range(L.NA) if av[a]]
        did = False
        for a in nz:
            for w in sorted(set(L.avars[a])):
                if w in FORBID: continue
                tgt = T.solve_lin(a, w, v)
                if tgt is None or tgt == v[w]: continue
                cands = []
                if w in FREE and len(L.var_atoms[w]) == 1:
                    cands.append((w, tgt))
                elif w not in FREE:
                    d = definer.get(w)
                    if d is not None:
                        vv = list(v); vv[w] = tgt
                        for u in sorted(set(L.avars[d])):
                            if (u != w and u in FREE and u not in FORBID
                                    and len(L.var_atoms[u]) == 1):
                                nv = T.solve_lin(d, u, vv)
                                if nv is not None: cands.append((u, nv))
                for u, nv in cands:
                    v[u] = nv; did = True; break
                if did: break
            if did: break
        if not did: break
        ad.fwd(v, rounds=6)
    return v

v0 = L.load(os.path.join(HERE, 'mod9118_0.json'))
vm0 = [x % P for x in v0]
av0 = L.all_atom_values(v0)
g = {a: ad.grad(a, vm0) for a in BAD}
U = sorted((set(g[BAD[0]]) | set(g[BAD[1]])) - FORBID,
           key=lambda u: len(L.var_atoms[u]))
print(f'{len(U)} candidate free inputs; base score {score(v0)}', flush=True)
CAND = U[:34]
print(f'pairs over the {len(CAND)} lowest-consumer inputs: '
      f'{len(CAND)*(len(CAND)-1)//2}', flush=True)
r = [(-av0[a]) % P for a in BAD]
best = (score(v0), None)
t0 = time.time()
tested = solved = 0
for u1, u2 in itertools.combinations(CAND, 2):
    a11, a21 = g[BAD[0]].get(u1, 0) % P, g[BAD[1]].get(u1, 0) % P
    a12, a22 = g[BAD[0]].get(u2, 0) % P, g[BAD[1]].get(u2, 0) % P
    det = (a11 * a22 - a12 * a21) % P
    if det == 0: continue
    di = pow(det, -1, P)
    d1 = (a22 * r[0] - a12 * r[1]) % P * di % P
    d2 = (a11 * r[1] - a21 * r[0]) % P * di % P
    v = list(v0); v[u1] = v[u1] + d1; v[u2] = v[u2] + d2
    ad.fwd(v, rounds=6)
    av = L.all_atom_values(v)
    tested += 1
    if av[BAD[0]] % P or av[BAD[1]] % P:
        continue                      # linear model failed here
    solved += 1
    settle(v)
    s = score(v)
    if s > best[0]:
        best = (s, (u1, u2))
        T.save(v, os.path.join(HERE, f'pair_{s}.json'))
        av = L.all_atom_values(v)
        print(f'  *** x_{u1}+x_{u2}: score {s}  '
              f'nonzero {[a for a in range(L.NA) if av[a]]}', flush=True)
    if tested % 40 == 0:
        print(f'  {tested} tested, {solved} congruences solved, best {best[0]} '
              f'({time.time()-t0:.0f}s)', flush=True)
print(f'\ntested {tested}, both congruences solved in {solved}')
print(f'BEST {best[0]} via {best[1]}')
