"""S11 step 64: the CONSTANT-ROUTING network, as a unification problem.

Every blocker in this instance has one of two shapes:
    c*(A - B) - wire*handle        =>  A == B  (mod p)      [equality gadget]
    b*(x - K) - c*z   /   x - K    =>  x == K  (mod p)      [constant pin, gated by b]
Values are routed through a MUX, so which constant reaches which assertion depends
on the boolean controls.  Union-find the UNCONDITIONAL assertions first: if they
alone force x == K1 and x == K2 with K1 != K2, the instance is infeasible outright,
no branch can save it.  If not, the propagation says which controls are needed.
"""
import os, sys, collections
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, ad
P = ad.P
v = L.load(os.path.join(LAB, 'best', 'new_instance_partial_39026.json'))
WIRE = {u for u in range(L.NVARS) if v[u] == P}
print(f'wire members (value = p): {len(WIRE)}')

eqs, pins, cond = [], [], []
for a in range(L.NA):
    poly = L.polys[a]
    lins = [(m[0], c) for m, c in poly.items() if len(m) == 1]
    quads = [(m, c) for m, c in poly.items() if len(m) == 2]
    const = poly.get((), 0)
    # equality gadget: two opposite linear terms + a quadratic with a wire factor
    if len(lins) == 2 and len(quads) <= 1 and not const:
        (u1, c1), (u2, c2) = lins
        if c1 == -c2 and quads:
            (qm, qc) = quads[0]
            if qm[0] in WIRE or qm[1] in WIRE:
                eqs.append((u1, u2, a)); continue
    if len(lins) == 3 and not quads and not const:
        s = sorted(lins, key=lambda t: -abs(t[1]))
        if s[0][1] == -s[1][1] and abs(s[2][1]) == 1:
            eqs.append((s[0][0], s[1][0], a)); continue
    # unconditional constant pin: x - K
    if len(lins) == 1 and not quads and const and abs(const) > 2**40:
        pins.append((lins[0][0], (-const * pow(lins[0][1], -1, P)) % P, a)); continue
    # gated pin: b*(x - K) - c*z  ->  the b-linear coefficient is -K
    if quads and len(lins) >= 1 and not const:
        for (qm, qc) in quads:
            for b in qm:
                other = qm[0] if qm[1] == b else qm[1]
                kc = poly.get((b,), 0)
                if kc and qc and other != b:
                    cond.append((other, (-kc * pow(qc, -1, P)) % P, b, a))
                    break
            break
print(f'equality gadgets: {len(eqs)}   unconditional pins: {len(pins)}   '
      f'gated pins: {len(cond)}')

par = {}
def find(x):
    par.setdefault(x, x)
    while par[x] != x:
        par[x] = par[par[x]]; x = par[x]
    return x
def union(a, b):
    ra, rb = find(a), find(b)
    if ra != rb: par[ra] = rb
for u1, u2, a in eqs: union(u1, u2)
print(f'\nunion-find over the {len(eqs)} equality gadgets: '
      f'{len(set(find(x) for x in par))} classes over {len(par)} variables')

anchor = {}
clash = []
for u, K, a in pins:
    r = find(u)
    if r in anchor and anchor[r][0] != K: clash.append((u, K, anchor[r], a))
    else: anchor[r] = (K, u, a)
print(f'unconditional pins anchor {len(anchor)} classes; CLASHES: {len(clash)}')
for c in clash[:6]: print(f'   x_{c[0]} pinned to {str(c[1])[:24]}... but its class '
                          f'already anchored to {str(c[2][0])[:24]}... (atoms {c[3]}, {c[2][2]})')
if clash:
    print('\n*** UNCONDITIONAL CLASH -> the instance would be infeasible outright')
else:
    print('\nno unconditional clash: every equality class carries at most one constant')
    print('   => the obstruction must involve the GATED pins, i.e. the branch choice')
# how many gated pins land on an already-anchored class with a different constant?
bad = 0
for u, K, b, a in cond:
    r = find(u)
    if r in anchor and anchor[r][0] != K: bad += 1
print(f'gated pins whose constant disagrees with their class anchor: {bad} '
      f'of {len(cond)}')
