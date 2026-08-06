"""S10 step 52: apply wire-deformation kernel vectors and measure.

The kernel guarantees only that the wire-IDENTITY equations stay satisfied.
Everything else (products wire*handle and their downstream checks) must be
re-verified empirically.  Also: does the kernel move the root pin x_26064?
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
print(f'wire {len(WIRE)} members, {len(BASIS)} kernel directions')

ROOT = 26064
print(f'\nroot x_{ROOT} in wire: {ROOT in widx}')
if ROOT in widx:
    j = widx[ROOT]
    print('  kernel coefficients at the root:')
    for i, b in enumerate(BASIS):
        print(f'    basis {i}: {str(b[j])[:60]}... ({len(str(abs(b[j])))} digits)' if b[j]
              else f'    basis {i}: 0')
    g = 0
    for b in BASIS:
        g = math.gcd(g, abs(b[j]))
    print(f'  gcd at root = {"0 (root is FIXED by the kernel)" if g == 0 else str(g)[:40]+f" ({len(str(g))} digits)"}')

# which members are moved at all?
moved = [WIRE[j] for j in range(len(WIRE)) if any(b[j] for b in BASIS)]
print(f'\nmembers the kernel can move: {len(moved)} of {len(WIRE)}')
for u in (15616, 11360, 26064, 28599, 17499, 22665, 28961):
    j = widx.get(u)
    if j is None:
        print(f'  x_{u}: not a wire member')
        continue
    g = 0
    for b in BASIS:
        g = math.gcd(g, abs(b[j]))
    print(f'  x_{u:<7} gcd={"0 -> FIXED" if g == 0 else str(len(str(g)))+" digits"}   '
          f'{"(handle multiplier)" if u in (15616,11360,28599,17499,22665,28961) else ""}')

base = L.load(os.path.join(HERE, 'forward_state.json'))
av0 = L.all_atom_values(base)
f0 = L.failing_eqs(av0)
print(f'\nbase: failing={len(f0)} score={L.NEQ-len(f0)}')

for i, b in enumerate(BASIS):
    for scale in (1, -1):
        v = list(base)
        for j, u in enumerate(WIRE):
            v[u] = v[u] + scale * b[j]
        # wire members are gate outputs; block their definers so the deformation holds
        block = set()
        for u in WIRE:
            dd = L.definer.get(u)
            if dd is not None:
                block.add(dd)
        for _ in range(4):
            for x in ad.ORDER:
                a = L.definer[x]
                if a in block:
                    continue
                nv = T.solve_lin(a, x, v)
                if nv is not None:
                    v[x] = nv
        av = L.all_atom_values(v)
        nz = [a for a in range(L.NA) if av[a]]
        fail = L.failing_eqs(av)
        print(f'  basis {i} scale {scale:+d}: nonzero atoms={len(nz)} failing={len(fail)} '
              f'score={L.NEQ-len(fail)}')
