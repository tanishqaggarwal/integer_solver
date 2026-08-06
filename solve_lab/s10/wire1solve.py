"""S10 step 43: solve on the wire=1 branch, where every handle is unquantised.

With the wire at 1 the handles have granularity 1, so each of the congruence
checks becomes an exact integer equation with a free solo handle to absorb it:

  a7930  = 9367949*(x_24548 - x_25442) - x_7927 ,  x_7927  = wire*x_11052
  a29539 = 12846437*(x_14853 - x_1308) - x_29967,  x_29967 = wire*x_30163
  a35759 = 5113045*x_9118 - x_29854             ,  x_29854 = wire*x_1329
  a35760 = x_31864 - wire*x_10903

Close everything reachable, then report exactly what survives.
"""
import os, sys, json, collections
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L
import tools as T
import ad

P = ad.P
atom_out = L.atom_out
ROOTATOM = 37694
BLOCK = {ROOTATOM}

# solo free handles: free inputs occurring in exactly one atom
SOLO = {}
for u in range(L.NVARS):
    if u in L.definer:
        continue
    ats = L.var_atoms[u]
    if len(ats) == 1:
        SOLO.setdefault(ats[0], []).append(u)


def fwd_block(v, rounds=4):
    for _ in range(rounds):
        for u in ad.ORDER:
            a = L.definer[u]
            if a in BLOCK:
                continue
            nv = T.solve_lin(a, u, v)
            if nv is not None:
                v[u] = nv
    return v


def status(v):
    av = L.all_atom_values(v)
    nz = [a for a in range(L.NA) if av[a]]
    fail = L.failing_eqs(av)
    return av, nz, fail


v = L.load(os.path.join(HERE, 'wire1_state.json'))
av, nz, fail = status(v)
print(f'start: nz={nz} failing={len(fail)} score={L.NEQ-len(fail)}', flush=True)

for it in range(30):
    av, nz, fail = status(v)
    todo = [a for a in nz if a != ROOTATOM]
    print(f'\niter {it}: nz={nz} failing={len(fail)} score={L.NEQ-len(fail)}', flush=True)
    if not todo:
        print('  only the wire pin a37694 remains', flush=True)
        break
    progressed = False
    for a in todo:
        # candidates: solo handles anywhere in this atom's defining chain, then any var
        cands = []
        for u in sorted(L.avars[a]):
            if u in L.definer:
                d = L.definer[u]
                for h in SOLO.get(d, []):
                    cands.append((0, h, d))
            if u in SOLO.get(a, []):
                cands.append((0, u, a))
        for u in sorted(L.avars[a]):
            cands.append((1, u, a))
        for pri, u, tgt in cands:
            w = list(v)
            if tgt == a:
                nv = T.solve_lin(a, u, w)
            else:
                # solve the chain: pick handle value making atom a vanish
                lo = T.lin_parts(a, L.atom_out[L.definer[u]][1] if False else u, w)
                nv = None
                c0, r0 = (T.lin_parts(a, u, w) or (0, 0))
                if c0:
                    if r0 % c0 == 0:
                        nv = -r0 // c0
            if nv is None or nv == w[u]:
                continue
            w[u] = nv
            fwd_block(w)
            av2, nz2, fail2 = status(w)
            if len(nz2) < len(nz) or (len(nz2) == len(nz) and len(fail2) < len(fail)):
                print(f'   a{a}: set x_{u} -> nz={nz2} failing={len(fail2)}', flush=True)
                v = w; progressed = True
                break
        if progressed:
            break
    if not progressed:
        print('  no improving handle move', flush=True)
        break

av, nz, fail = status(v)
print(f'\nFINAL nz={nz} failing={len(fail)} score={L.NEQ-len(fail)}')
print(f'failing eqs: {fail[:30]}')
T.save(v, os.path.join(HERE, 'wire1_solved.json'))
