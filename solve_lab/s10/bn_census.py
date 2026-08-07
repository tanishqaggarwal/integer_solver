"""bn_census: enumerate boolean atoms x*x - x, classify their variables."""
import os, sys, json, collections
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad

P = 2**256 - 2**32 - 977

print('NA', L.NA, 'NEQ', L.NEQ, 'NVARS', L.NVARS)
print('FREE inputs', len(ad.FREE))

# --- find boolean atoms -------------------------------------------------
# poly is dict {monomial tuple: coeff}.  boolean = {(u,u):c, (u,):-c}
bools = {}          # atom -> (var, coeff)  meaning poly = c*(x^2 - x)
other_quad = []
for a in range(L.NA):
    Pp = L.polys[a]
    if len(Pp) != 2:
        continue
    ks = sorted(Pp.keys(), key=len)
    if len(ks[0]) == 1 and len(ks[1]) == 2 and ks[1][0] == ks[1][1] == ks[0][0]:
        u = ks[0][0]
        c1 = Pp[ks[0]]; c2 = Pp[ks[1]]
        if c1 == -c2:
            bools[a] = (u, c2)
        else:
            other_quad.append((a, u, c1, c2))

print('boolean atoms (c*(x^2-x)):', len(bools))
print('near-boolean with mismatched coeffs:', len(other_quad), other_quad[:10])
bvars = set(u for u,_ in bools.values())
print('distinct boolean vars:', len(bvars))
cc = collections.Counter(c for _,c in bools.values())
print('coefficient histogram:', cc.most_common(10))

FREESET = set(ad.FREE)
freeb = [u for u in bvars if u in FREESET]
print('boolean vars that are FREE inputs:', len(freeb))

# --- classify each boolean var ------------------------------------------
var2bool = collections.defaultdict(list)
for a,(u,c) in bools.items():
    var2bool[u].append(a)

rows = []
for u in sorted(bvars):
    atoms = L.var_atoms[u]           # all atoms mentioning u
    ba = var2bool[u]
    others = [a for a in atoms if a not in bools]
    # equations touched by the boolean atoms
    beq = set()
    for a in ba: beq.update(L.atom2eq.get(a, ()))
    # equations touched by other atoms of u
    oeq = set()
    for a in others: oeq.update(L.atom2eq.get(a, ()))
    rows.append(dict(u=u, free=(u in FREESET), nbool=len(ba), nother=len(others),
                     beq=len(beq), oeq=len(oeq), bool_atoms=ba,
                     defined=(u in L.definer)))

clean = [r for r in rows if r['nother']==0]
cleanfree = [r for r in clean if r['free']]
print()
print('boolean vars appearing in NO other atom (clean carriers):', len(clean))
print('   of which FREE inputs:', len(cleanfree))
h = collections.Counter(r['nother'] for r in rows)
print('nother histogram:', sorted(h.items())[:15])
h2 = collections.Counter((r['free'], r['nother']==0) for r in rows)
print('(free, clean) histogram:', h2)

json.dump({'bools': {str(a): list(v) for a,v in bools.items()},
           'rows': rows}, open(os.path.join(HERE,'bn_census.json'),'w'))
print('saved bn_census.json')
