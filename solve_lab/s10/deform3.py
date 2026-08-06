"""S10 step 55: adversarial re-test of my own section-18 closure.

Section 18 concluded the wire is pinned to +-p because every member appears in a
monomial w_i*w_j and invariance of that MONOMIAL forces w_i*w_j = p^2.  That is a
SUFFICIENT condition, not a necessary one: the square check is E^2 = 0 <=> E = 0,
and E contains other variables -- possibly FREE handles that can absorb a change
in w_u.  So: deform along the kernel, restore the product gates, and then try to
repair every remaining broken atom through the free variables inside it.
"""
import os, sys, json, collections, time
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
_d = json.load(open(os.path.join(HERE, 'wirekernel.json')))
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L
import tools as T
import ad

P = ad.P
WIRE, BASIS = _d['wire'], _d['basis']
WSET = set(WIRE)
base = L.load(os.path.join(HERE, 'forward_state.json'))
WIREDEF = set(L.definer[u] for u in WIRE if u in L.definer)

PROD = []
for a, out in L.atom_out.items():
    t = out[1]; vs = L.avars[a]
    wm = [u for u in vs if u in WSET]
    fr = [u for u in vs if u not in L.definer]
    if len(vs) == 3 and len(wm) == 1 and len(fr) == 1 and t != wm[0] and t != fr[0]:
        PROD.append((a, t, wm[0], fr[0]))
PRODA = set(a for a, t, w, f in PROD)
BLOCK = WIREDEF | PRODA


def fwd(v, rounds=3):
    for _ in range(rounds):
        for x in ad.ORDER:
            a = L.definer[x]
            if a in BLOCK:
                continue
            nv = T.solve_lin(a, x, v)
            if nv is not None:
                v[x] = nv
    return v


def build(coeffs):
    v = list(base)
    for j, u in enumerate(WIRE):
        v[u] += sum(c * b[j] for c, b in zip(coeffs, BASIS))
    for a, t, wm, fr in PROD:
        w = v[wm]
        if w and base[t] % w == 0:
            v[fr] = base[t] // w
            v[t] = base[t]
    fwd(v)
    return v


f0 = len(L.failing_eqs(L.all_atom_values(base)))
print(f'base failing={f0} score={L.NEQ-f0}', flush=True)

for coeffs in ([1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 0], [1, -1, 0], [1, 1, 1]):
    v = build(coeffs)
    av = L.all_atom_values(v)
    nz = [a for a in range(L.NA) if av[a]]
    hard = [a for a in nz if a not in WIREDEF]
    fail = L.failing_eqs(av)
    print(f'\ncoeffs {coeffs}: nz={len(nz)} (non-copy {len(hard)}) '
          f'failing={len(fail)} score={L.NEQ-len(fail)}', flush=True)

    # try to repair each non-copy broken atom through a free variable inside it
    repaired = []
    for _ in range(len(hard) + 4):
        av = L.all_atom_values(v)
        hard = [a for a in range(L.NA) if av[a] and a not in WIREDEF]
        if not hard:
            break
        best = None
        for a in hard:
            for u in sorted(L.avars[a]):
                if u in L.definer or u in WSET:
                    continue
                nv = T.solve_lin(a, u, v)
                if nv is None or nv == v[u]:
                    continue
                w = list(v); w[u] = nv; fwd(w)
                aw = L.all_atom_values(w)
                nh = len([x for x in range(L.NA) if aw[x] and x not in WIREDEF])
                nf = len(L.failing_eqs(aw))
                if best is None or (nh, nf) < best[0]:
                    best = ((nh, nf), a, u, w)
        if best is None or best[0][0] >= len(hard):
            break
        (nh, nf), a, u, w = best
        repaired.append((a, u))
        v = w
    av = L.all_atom_values(v)
    hard = [a for a in range(L.NA) if av[a] and a not in WIREDEF]
    fail = L.failing_eqs(av)
    print(f'   after repair ({len(repaired)} moves): non-copy broken={len(hard)} '
          f'{hard[:14]} failing={len(fail)} score={L.NEQ-len(fail)}', flush=True)
    if len(fail) < f0:
        T.save(v, os.path.join(HERE, f'deform3_{"".join(map(str,coeffs))}.json'))
