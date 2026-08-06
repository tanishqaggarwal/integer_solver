"""S10 step 80: slide along the cyclic components' kernels.

40 two-variable SCCs are under-determined (kernel dim 1): each admits a
one-parameter family of solutions with BOTH its gate atoms still satisfied.
Forward evaluation silently picks one point on each line; no local analysis of
mine could see the rest.  That is 40 genuinely new free parameters.

For each, slide along the line and measure: does anything change, and in
particular do the binding residues D0 / K2 move?  A slide that moves either
residue kills a congruence and beats the deliverable.
"""
import os, sys, collections, json
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L
import tools as T

P = 2**256 - 2**32 - 977
spec = json.load(open(os.path.join(HERE, 'cycles.json')))
UNDER = spec['under']
print(f'under-determined components: {len(UNDER)}')

base = L.load(os.path.join(LAB, 'best', 'new_instance_partial_39026.json'))
bav = L.all_atom_values(base)
BASE_NZ = set(a for a in range(L.NA) if bav[a])
f0 = L.failing_eqs(bav)
D0 = (base[7068] - base[2099]) % P
K2 = base[28730] % P
print(f'witness failing {len(f0)}; D0 mod p = {str(D0)[:24]}...; '
      f'K2 = {str(K2)[:24]}...')
BLOCK = set(BASE_NZ) | {22231}

useful = []
for kdim, msize, cs in UNDER:
    a, b = cs[0], cs[1]
    ats = sorted(set(L.var_atoms[a]) | set(L.var_atoms[b]))
    gates = [L.definer.get(a), L.definer.get(b)]
    extra = [x for x in ats if x not in gates]
    # slide: set x_a = v[a] + d, then re-solve x_b from a's OWN gate partner
    for d in (1, 2, P):
        w = list(base)
        w[a] = base[a] + d
        # re-solve x_b from its gate atom, then re-solve x_a from its own,
        # iterating: on a degenerate pair this converges to a NEW valid point
        ok = True
        for _ in range(6):
            nb = T.solve_lin(L.definer[b], b, w)
            if nb is None: ok = False; break
            w[b] = nb
        if not ok:
            continue
        # verify both gate atoms still vanish
        ga = L.evalpoly(L.polys[L.definer[a]], w) if hasattr(L, 'evalpoly') else None
        gb = L.evalpoly(L.polys[L.definer[b]], w) if hasattr(L, 'evalpoly') else None
        if ga or gb:
            continue
        try:
            L.ripple(w, {}, block=BLOCK)
        except Exception:
            pass
        aw = L.all_atom_values(w)
        nz = set(x for x in range(L.NA) if aw[x])
        newbad = sorted(nz - BASE_NZ)
        fail = L.failing_eqs(aw)
        nD = (w[7068] - w[2099]) % P
        nK = w[28730] % P
        moved = (nD != D0) or (nK != K2)
        if not newbad or moved:
            useful.append((len(fail), a, b, d, newbad, moved))
            print(f'  SCC ({a},{b}) d={"p" if d == P else d}: extra atoms {newbad}, '
                  f'failing {len(fail)}, residues moved: {moved}', flush=True)
        break

print(f'\nslides with no collateral or a moved residue: {len(useful)}')
if useful:
    useful.sort()
    for f, a, b, d, nb, mv in useful[:12]:
        print(f'  ({a},{b}) failing={f} score={L.NEQ-f} newbad={nb} residue_moved={mv}')
else:
    print('  none -- every cyclic slide either breaks atoms or is invisible')

print('\n=== how many cyclic pairs touch ONLY their own two gate atoms? ===')
inert = 0
for kdim, msize, cs in UNDER:
    a, b = cs[0], cs[1]
    ats = set(L.var_atoms[a]) | set(L.var_atoms[b])
    if ats <= {L.definer.get(a), L.definer.get(b)}:
        inert += 1
print(f'  {inert} of {len(UNDER)} are inert (sliding changes literally nothing)')
