"""S10 step 63: price every wire member as a deformation target.

Two known routes and their costs:
   kernel deformation : 0 identity cost, support 215 -> ~20 square checks break
   uniform wire = 1   : squares fine, but the root pin a37694 costs 12 (+1) = 13

Intermediate route: move ONE member (plus its copy-descendants).  Cost =
   price(copy atom that defines it)  +  |square-check equations containing it|
Budget to beat: 7.  Report the cheapest members, and flag those that are handle
multipliers for certificate 1's checks (x_14466, x_3915, x_11360).
"""
import os, sys, collections, json
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L
import tools as T
import ad

P = ad.P
atom_out = L.atom_out
base = L.load(os.path.join(HERE, 'forward_state.json'))
WIRE = sorted(u for u in range(L.NVARS) if base[u] == P)
WSET = set(WIRE)

# square/high-degree checks containing each member with multiplicity >= 2
sqeq = collections.defaultdict(set)
for a in range(L.NA):
    if a in atom_out:
        continue
    for m, c in L.polys[a].items():
        cnt = collections.Counter(u for u in m if u in WSET)
        for u, k in cnt.items():
            if k >= 2:
                sqeq[u] |= set(L.atom2eq.get(a, ()))

# copy-descendants of a member (members defined FROM it)
def descend(u):
    out, fr = {u}, [u]
    while fr:
        nxt = []
        for w in fr:
            for a in L.var_atoms[w]:
                if a not in atom_out:
                    continue
                t = atom_out[a][1]
                if t != w and t in WSET and t not in out and all(z in WSET for z in L.avars[a]):
                    out.add(t); nxt.append(t)
        fr = nxt
    return out

HANDLE_MULT = {14466, 3915, 11360, 15616, 28599, 17499, 22665, 28961}
rows = []
for u in WIRE:
    dset = descend(u)
    ident = set()
    for t in dset:
        d = L.definer.get(t)
        if d is not None:
            ident |= set(L.atom2eq.get(d, ()))
    sq = set()
    for t in dset:
        sq |= sqeq.get(t, set())
    total = len(ident | sq)
    rows.append((total, len(ident), len(sq), len(dset), u))
rows.sort()
print(f'{"total":>6} {"ident":>6} {"square":>7} {"desc":>5}  member')
for tot, ni, ns, nd, u in rows[:25]:
    tag = '  <== HANDLE MULTIPLIER' if u in HANDLE_MULT else ''
    print(f'{tot:>6} {ni:>6} {ns:>7} {nd:>5}  x_{u}{tag}')

print('\ncost for the certificate-1 handle multipliers specifically:')
for u in sorted(HANDLE_MULT):
    for tot, ni, ns, nd, w in rows:
        if w == u:
            print(f'  x_{u:<7} total={tot:<5} identity={ni:<4} square={ns:<4} '
                  f'descendants={nd}')
            break

print(f'\ncheapest member overall: {rows[0][4]} at {rows[0][0]} equations '
      f'(budget to beat: 7)')
json.dump([{'member': u, 'total': t, 'ident': i, 'square': s, 'desc': d}
           for t, i, s, d, u in rows],
          open(os.path.join(HERE, 'memberprice.json'), 'w'))
