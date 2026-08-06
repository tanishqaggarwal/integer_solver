"""S10 step 30: what ACTUALLY imposes the two congruences?

Under the canonical orientation, five of the seven residual atoms are GATES:
  a22229 -> x_7068,  a22230 -> x_28730,  a35758 -> x_29854,
  a35761 -> x_31864, a35762 -> x_642
so they vanish by construction.  Only a35759 and a35760 are checks there.
Therefore the two congruences are NOT intrinsic to these atoms -- they are imposed
by DOWNSTREAM checks.  Identify them exactly.
"""
import os, sys, collections
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L
import tools as T

P = 2**256 - 2**32 - 977
NZ = [22229, 22230, 35758, 35759, 35760, 35761, 35762]
v = L.load(os.path.join(LAB, 'best', 'new_instance_partial_39026.json'))

print('=== orientation of the seven residual atoms ===')
for a in NZ:
    out = L.atom_out.get(a)
    kind = f'GATE -> defines x_{out[1]} (coeff {out[0]})' if out else 'CHECK'
    print(f'  a{a}: {kind}')
    print(f'      {L.atom_src[a]}')

print('\n=== the two downstream checks that actually bind ===')
for a in (29539, 7930):
    print(f'\na{a}: {L.atom_src[a]}')
    out = L.atom_out.get(a)
    print(f'   orientation: {"GATE -> x_%d" % out[1] if out else "CHECK"}')
    print(f'   in {len(L.atom2eq.get(a,{}))} equations')
    for u in sorted(L.avars[a]):
        d = L.definer.get(u)
        print(f'     x_{u:<7} {"FREE" if d is None else "def by a%d" % d:<16} '
              f'natoms={len(L.var_atoms[u]):<3} val%p={v[u] % P if v[u] else 0}')
        if d is not None:
            print(f'         {L.atom_src[d][:120]}')

# Trace: which free inputs does x_1308 depend on?  and x_25442?
print('\n=== dependency cone (free inputs only) of the binding wires ===')
avars = L.avars


def cone(target, maxdepth=40):
    """Free inputs reachable backwards from `target` through gate definitions."""
    seen, free, frontier = {target}, set(), [target]
    for _ in range(maxdepth):
        nxt = []
        for u in frontier:
            d = L.definer.get(u)
            if d is None:
                free.add(u); continue
            for w in avars[d]:
                if w != u and w not in seen:
                    seen.add(w); nxt.append(w)
        if not nxt:
            break
        frontier = nxt
    return free, seen


for t in (1308, 25442, 12186, 22649, 2099, 24548, 14853):
    f, s = cone(t)
    print(f'  x_{t:<7} cone: {len(s)} vars, {len(f)} free inputs')
