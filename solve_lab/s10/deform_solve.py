"""S10 step 62: deform the wire, then USE the freed handles to solve the checks.

deform2/deform3 restored each handle so the product gate reproduced its ORIGINAL
output -- that preserves the state but wastes the whole point of the deformation.
The right move: after deforming, every handle has granularity 1, so SOLVE the
failing checks through them.

Certificate 1's members all carry a p-wire multiplier in their handle:
    a2423   x_9899  = x_14466 * x_14768
    a31670  x_29309 = x_105   * x_3915
    a29539  x_29967 = x_11360 * x_30163
so the wire deformation is precisely the lever on the one expensive certificate.
"""
import os, sys, json, collections, math, time
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
_k = json.load(open(os.path.join(HERE, 'wirekernel.json')))
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L
import tools as T
import ad

P = ad.P
atom_out = L.atom_out
WIRE, BASIS = _k['wire'], _k['basis']
WSET = set(WIRE)
base = L.load(os.path.join(HERE, 'forward_state.json'))
WIREDEF = set(L.definer[u] for u in WIRE if u in L.definer)

SOLO = collections.defaultdict(list)
for u in range(L.NVARS):
    if u not in L.definer and len(L.var_atoms[u]) == 1:
        SOLO[L.var_atoms[u][0]].append(u)


def fwd(v, rounds=3):
    for _ in range(rounds):
        for x in ad.ORDER:
            a = L.definer[x]
            if a in WIREDEF:
                continue
            nv = T.solve_lin(a, x, v)
            if nv is not None:
                v[x] = nv
    return v


def status(v):
    av = L.all_atom_values(v)
    nz = [a for a in range(L.NA) if av[a]]
    return av, nz, L.failing_eqs(av)


def handle_moves(v, a):
    """Ways to zero atom a through a solo free handle one gate upstream."""
    out = []
    for t in sorted(L.avars[a]):
        r = T.lin_parts(a, t, v)
        if r is None:
            continue
        c, rest = r
        if c == 0 or rest % c:
            continue
        tn = -rest // c
        if tn == v[t]:
            continue
        if t not in L.definer:
            out.append((t, tn)); continue
        d = L.definer[t]
        w = list(v); w[t] = tn
        for h in SOLO.get(d, []):
            nv = T.solve_lin(d, h, w)
            if nv is not None and nv != v[h]:
                out.append((h, nv))
    return out


f0 = len(L.failing_eqs(L.all_atom_values(base)))
print(f'base failing={f0} score={L.NEQ-f0}', flush=True)
t0 = time.time()

for coeffs in ([1, 0, 0], [0, 1, 0], [0, 0, 1]):
    v = list(base)
    for j, u in enumerate(WIRE):
        v[u] += sum(c * b[j] for c, b in zip(coeffs, BASIS))
    fwd(v)
    av, nz, fail = status(v)
    print(f'\ncoeffs {coeffs}: after deform nz={len(nz)} failing={len(fail)} '
          f'score={L.NEQ-len(fail)} ({time.time()-t0:.0f}s)', flush=True)
    print(f'   granularity check: |w_3915|={len(str(abs(v[3915])))}d '
          f'|w_11360|={len(str(abs(v[11360])))}d', flush=True)

    # now SOLVE every nonzero non-copy atom through its (freed) handle
    for rnd in range(40):
        av, nz, fail = status(v)
        todo = [a for a in nz if a not in WIREDEF]
        if not todo:
            break
        best = None
        for a in todo:
            for h, val in handle_moves(v, a):
                w = list(v); w[h] = val; fwd(w)
                aw, nzw, fw = status(w)
                todow = len([x for x in nzw if x not in WIREDEF])
                key = (len(fw), todow)
                if best is None or key < best[0]:
                    best = (key, a, h, w)
        if best is None or best[0][0] >= len(fail):
            break
        key, a, h, w = best
        v = w
        print(f'   rnd {rnd}: closed a{a} via handle x_{h} -> failing={key[0]} '
              f'score={L.NEQ-key[0]}', flush=True)
    av, nz, fail = status(v)
    todo = [a for a in nz if a not in WIREDEF]
    print(f'   FINAL non-copy broken={len(todo)} {todo[:14]} failing={len(fail)} '
          f'score={L.NEQ-len(fail)}', flush=True)
    if len(fail) < f0:
        T.save(v, os.path.join(HERE, f'defsolve_{L.NEQ-len(fail)}.json'))
        print(f'   saved defsolve_{L.NEQ-len(fail)}.json', flush=True)
