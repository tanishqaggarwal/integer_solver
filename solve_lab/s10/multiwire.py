"""S10 step 54: the 13 multi-wire monomials — the last obstruction on the crack.

Their invariance condition  p(d_i + d_j) + d_i*d_j = 0  factors as

        w_i * w_j = p^2

so it does NOT force w = p: setting w_i = 1 and w_j = p^2 satisfies it.  Find the
actual pairs, see whether the constraint graph is consistent, and check whether it
is compatible with the 3-dimensional linear kernel.
"""
import os, sys, json, collections, math
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

pairs = collections.Counter()
atoms = []
for a in range(L.NA):
    for m, c in L.polys[a].items():
        wm = [u for u in m if u in WSET]
        if len(wm) >= 2:
            pairs[tuple(sorted(wm[:2]))] += 1
            atoms.append((a, tuple(sorted(wm[:2])), m, c))
print(f'atoms containing a multi-wire monomial: {len(set(x[0] for x in atoms))}')
print(f'distinct wire pairs: {len(pairs)}')
for pr, n in pairs.most_common():
    i, j = pr
    gi = math.gcd(*[abs(b[widx[i]]) for b in BASIS]) if len(BASIS) > 1 else abs(BASIS[0][widx[i]])
    gj = math.gcd(*[abs(b[widx[j]]) for b in BASIS]) if len(BASIS) > 1 else abs(BASIS[0][widx[j]])
    gi = 0
    gj = 0
    for b in BASIS:
        gi = math.gcd(gi, abs(b[widx[i]])); gj = math.gcd(gj, abs(b[widx[j]]))
    print(f'  (x_{i}, x_{j})  in {n} monomials   gcd_i={gi} gcd_j={gj}')

print('\nthe atoms involved:')
seen = set()
for a, pr, m, c in atoms:
    if a in seen:
        continue
    seen.add(a)
    out = L.atom_out.get(a)
    print(f'  a{a} [{"GATE->x_%d" % out[1] if out else "CHECK"}] neq={len(L.atom2eq.get(a,{}))}'
          f'  pair={pr}')
    print(f'      {L.atom_src[a][:170]}')

# constraint graph: w_i * w_j = p^2  on each edge
print('\n=== constraint graph  w_i * w_j = p^2 ===')
g = collections.defaultdict(set)
for (i, j) in pairs:
    g[i].add(j); g[j].add(i)
seen = set(); comps = []
for u in g:
    if u in seen: continue
    comp, st = set(), [u]
    while st:
        x = st.pop()
        if x in comp: continue
        comp.add(x); seen.add(x)
        st.extend(g[x])
    comps.append(comp)
print(f'components: {len(comps)}  sizes {[len(c) for c in comps]}')
for c in comps:
    # 2-colour: w alternates between t and p^2/t
    col = {}
    root = next(iter(c)); col[root] = 0; st = [root]; bipartite = True
    while st:
        x = st.pop()
        for y in g[x]:
            if y not in col:
                col[y] = 1 - col[x]; st.append(y)
            elif col[y] == col[x]:
                bipartite = False
    selfloop = any(i == j for (i, j) in pairs if i in c)
    print(f'  component {sorted(c)[:8]}... size {len(c)} bipartite={bipartite} '
          f'self-loop={selfloop}')
    if selfloop:
        print('     -> a self-pair w_i^2 = p^2 forces w_i = +-p')
