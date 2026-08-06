"""S10 step 53: deform the wire, then RE-SOLVE the handles.

Applying a kernel vector changes every product  x_B = wire * handle, so the gate
outputs move and downstream checks break.  But each handle is a free input, so we
can re-solve it to restore the ORIGINAL x_B -- and that is exactly divisible when
the new wire value divides the target.  Test the whole pipeline.
"""
import os, sys, json, math, collections
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L
import tools as T
import ad

P = ad.P
d = json.load(open(os.path.join(HERE, 'wirekernel.json')))
WIRE, BASIS = d['wire'], d['basis']
widx = {u: i for i, u in enumerate(WIRE)}
WSET = set(WIRE)
base = L.load(os.path.join(HERE, 'forward_state.json'))
base_av = L.all_atom_values(base)
WIREDEF = set(L.definer[u] for u in WIRE if u in L.definer)

# product gates of the form  x_B - wire * handle   (handle a free input)
PROD = []
for a, out in L.atom_out.items():
    t = out[1]
    vs = L.avars[a]
    wm = [u for u in vs if u in WSET]
    fr = [u for u in vs if u not in L.definer]
    if len(vs) == 3 and len(wm) == 1 and len(fr) == 1 and t != wm[0] and t != fr[0]:
        PROD.append((a, t, wm[0], fr[0]))
print(f'product gates  x_B = wire * free_handle : {len(PROD)}')


def apply_deform(coeffs):
    v = list(base)
    for j, u in enumerate(WIRE):
        dv = sum(c * b[j] for c, b in zip(coeffs, BASIS))
        v[u] = v[u] + dv
    return v


def resolve_handles(v):
    """Restore each product gate's ORIGINAL output by re-solving its handle."""
    fixed = failed = 0
    for a, t, wm, fr in PROD:
        target = base[t]
        w = v[wm]
        if w == 0:
            failed += 1; continue
        if target % w == 0:
            v[fr] = target // w
            v[t] = target
            fixed += 1
        else:
            failed += 1
    return fixed, failed


def fwd(v, rounds=3):
    block = WIREDEF | set(a for a, t, wm, fr in PROD)
    for _ in range(rounds):
        for x in ad.ORDER:
            a = L.definer[x]
            if a in block:
                continue
            nv = T.solve_lin(a, x, v)
            if nv is not None:
                v[x] = nv
    return v


f0 = L.failing_eqs(base_av)
print(f'base: failing={len(f0)} score={L.NEQ-len(f0)}')

for i in range(len(BASIS)):
    for s in (1, -1):
        co = [0] * len(BASIS); co[i] = s
        v = apply_deform(co)
        fx, fl = resolve_handles(v)
        fwd(v)
        av = L.all_atom_values(v)
        nz = [a for a in range(L.NA) if av[a]]
        fail = L.failing_eqs(av)
        print(f'  basis {i} x{s:+d}: handles restored {fx}, not divisible {fl}; '
              f'nz={len(nz)} failing={len(fail)} score={L.NEQ-len(fail)}', flush=True)
        # are the two hard handles now small?
        print(f'      w_15616 digits={len(str(abs(v[15616])))}  '
              f'w_11360 digits={len(str(abs(v[11360])))}')
