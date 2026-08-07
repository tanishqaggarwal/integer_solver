"""S11 step 65: unification WITH zero-pins.

x - handle*wire = 0  asserts  x == 0 (mod p).  That is the shape of the seven
residual atoms (a22230 pins x_28730 == 0, a35758 pins x_29854 == 0, ...), and the
first pass mislabelled it.  Add it, then propagate: a class forced to be both == 0
and == K != 0 is a contradiction, and the instance would be infeasible.
"""
import os, sys, collections
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, ad
P = ad.P
v = L.load(os.path.join(LAB, 'best', 'new_instance_partial_39026.json'))
WIRE = {u for u in range(L.NVARS) if v[u] == P}
BOOL = set()
for a, poly in enumerate(L.polys):
    ks = list(poly.items())
    if len(ks) == 2:
        sq = [m for m, c in ks if len(m) == 2 and m[0] == m[1]]
        li = [m for m, c in ks if len(m) == 1]
        if sq and li and sq[0][0] == li[0][0]: BOOL.add(li[0][0])

eqs, zeros, pins, gated = [], [], [], []
for a in range(L.NA):
    poly = L.polys[a]
    lins = [(m[0], c) for m, c in poly.items() if len(m) == 1]
    quads = [(m, c) for m, c in poly.items() if len(m) == 2 and m[0] != m[1]]
    const = poly.get((), 0)
    wq = [(m, c) for m, c in quads if m[0] in WIRE or m[1] in WIRE]
    if const == 0 and len(lins) == 1 and len(quads) == 1 and wq:
        zeros.append((lins[0][0], a)); continue                 # x == 0 (mod p)
    if const == 0 and len(lins) == 2 and len(quads) <= 1:
        (u1, c1), (u2, c2) = lins
        if c1 == -c2 and (not quads or wq):
            eqs.append((u1, u2, a)); continue                   # A == B (mod p)
    if const == 0 and len(lins) == 3 and not quads:
        s = sorted(lins, key=lambda t: -abs(t[1]))
        if s[0][1] == -s[1][1] and abs(s[2][1]) == 1:
            eqs.append((s[0][0], s[1][0], a)); continue
    if const and len(lins) == 1 and not quads and abs(const) > 2**40:
        pins.append((lins[0][0], (-const) * pow(lins[0][1], -1, P) % P, a)); continue
    if const == 0 and quads:
        for (qm, qc) in quads:
            b = qm[0] if qm[0] in BOOL else (qm[1] if qm[1] in BOOL else None)
            if b is None: continue
            other = qm[1] if qm[0] == b else qm[0]
            kc = poly.get((b,), 0)
            if kc:
                gated.append((other, (-kc) * pow(qc, -1, P) % P, b, a))
            break
print(f'equality gadgets {len(eqs)}   ZERO-pins {len(zeros)}   '
      f'constant pins {len(pins)}   gated pins {len(gated)}')

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
anchor, clash = {}, []
def put(u, K, a, kind):
    r = find(u)
    if r in anchor and anchor[r][0] != K:
        clash.append((u, K, kind, a, anchor[r]))
    else:
        anchor.setdefault(r, (K, u, a, kind))
for u, a in zeros: put(u, 0, a, 'zero')
for u, K, a in pins: put(u, K, a, 'const')
print(f'\nclasses {len(set(find(x) for x in par))}, anchored {len(anchor)}')
print(f'CLASHES among unconditional assertions: {len(clash)}')
for c in clash[:8]:
    print(f'   x_{c[0]} -> {str(c[1])[:22]} ({c[2]}, atom a{c[3]}) '
          f'but class anchored {str(c[4][0])[:22]} ({c[4][3]}, atom a{c[4][2]})')
if clash:
    print('\n*** UNCONDITIONAL CLASH: the instance is infeasible regardless of branch')
else:
    print('\nno unconditional clash')
g = 0
for u, K, b, a in gated:
    r = find(u)
    if r in anchor and anchor[r][0] != K: g += 1
print(f'gated pins disagreeing with their class anchor: {g} of {len(gated)}')
# and the key one: are the residual variables in anchored classes?
print('\nresidual variables and their class anchors:')
for u in (9118, 8731, 28730, 7068, 2099, 642, 29854, 31864):
    r = find(u)
    an = anchor.get(r)
    print(f'  x_{u:<7} class size {sum(1 for x in par if find(x) == r):<5} '
          f'anchor {("0" if an and an[0]==0 else str(an[0])[:20]) if an else "NONE"}'
          f'{" via " + an[3] + " a" + str(an[2]) if an else ""}')
