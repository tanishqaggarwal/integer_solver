"""S10 step 37: the certificate's atoms, and the CORE derived from scratch.

Certificate (s10/constrained.py): FAIL {7930,29539,35759,35760,40826,41512}
plus SAT {3576,3578,7938,7939,18691,18694} becomes inconsistent exactly when the
core check 19297 is added.

Derive the core's condition directly from the atoms instead of inheriting it,
then test the second branch (A*c^2 == B^2) as an explicit quadratic-residue
question, including how much freedom A actually has.
"""
import os, sys, collections
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L
import tools as T
import ad

P = ad.P
v = L.load(os.path.join(HERE, 'forward_state.json'))
vm = [x % P for x in v]
av = L.all_atom_values(v)

print('=== certificate atoms ===')
for a in (3576, 3578, 7938, 7939, 18691, 18694, 19297, 19299, 30984):
    out = L.atom_out.get(a)
    print(f'\na{a} [{"GATE->x_%d" % out[1] if out else "CHECK"}] '
          f'neq={len(L.atom2eq.get(a,{}))} value={str(av[a])[:40]}')
    print(f'   {L.atom_src[a][:300]}')

print('\n\n=== core variables, derived ===')
CORE = [29322, 3558, 33469, 1326, 27713, 15298, 14853, 12186, 24908, 16742,
        22162, 24453, 30213]
for u in CORE:
    d = L.definer.get(u)
    print(f'  x_{u:<7} {"FREE" if d is None else "def a%d" % d:<12} '
          f'natoms={len(L.var_atoms[u]):<3} val%p={v[u] % P}')
    if d is not None:
        print(f'      {L.atom_src[d][:120]}')

# --- the second branch as an explicit QR question --------------------------
print('\n=== second branch:  A*c^2 == B^2  (mod p) ? ===')
A = vm[33469]; c = vm[1326]; B = vm[27713]
u = vm[29322]; w = vm[3558]
print(f'  u = {u}\n  w = {w}')
print(f'  A = {A}\n  c = {c}\n  B = {B}')
print(f'  A*c^2 - B^2 mod p = {(A * c * c - B * B) % P}')
if c:
    tgt = B * pow(c, -1, P) % P
    tgt = tgt * tgt % P
    print(f'  required A = (B/c)^2 = {tgt}')
    print(f'  current  A = {A}')
    print(f'  delta needed on A = {(tgt - A) % P}')
leg = pow(A, (P - 1) // 2, P) if A else 0
print(f'  Legendre(A) = {1 if leg == 1 else (-1 if leg == P-1 else 0)}  '
      f'(A is a {"QR" if leg == 1 else "non-residue"})')

# --- how much freedom does A have? ----------------------------------------
print('\n=== sensitivity of A = x_33469 to free inputs (mod p) ===')
gA = ad.grad(33469, vm) if 33469 not in L.definer else None
if 33469 in L.definer:
    # A is a gate; use a pseudo-check: measure d(x_33469)/d(free) by AD on its definer
    print('  x_33469 is a gate; tracing its own cone instead')
    seen, frontier, freeins = {33469}, [33469], []
    for _ in range(40):
        nxt = []
        for t in frontier:
            d = L.definer.get(t)
            if d is None:
                freeins.append(t); continue
            for wv in L.avars[d]:
                if wv != t and wv not in seen:
                    seen.add(wv); nxt.append(wv)
        if not nxt: break
        frontier = nxt
    print(f'  A depends on {len(freeins)} free inputs: {[f"x_{x}" for x in freeins[:25]]}')
