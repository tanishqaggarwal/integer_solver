"""bn_struct: structure of the optimal deficiency block; boolean atoms in failing eqs."""
import os, sys, json, collections
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad

cen = json.load(open(os.path.join(HERE,'bn_census.json')))
bools = {int(a):tuple(t) for a,t in cen['bools'].items()}
BA = sorted(bools); FREESET=set(ad.FREE)
d = json.load(open(os.path.join(HERE,'bn_defic.json')))
S1 = d['all']['S']; E1 = set(d['all']['E'])
aeq = {a:set(L.atom2eq.get(a,{})) for a in BA}

FAIL = [12231, 12270, 12350, 14584, 18673, 22044, 29125]
print('=== boolean atoms inside the 7 failing equations ===')
for e in FAIL:
    m,sq,co = L.eq_atoms[e]
    bs = [(a,co[a]) for a in co if a in bools]
    print(f'  eq{e}: {len(co)} atoms, {len(bs)} boolean -> {bs}')
    for a,c in bs:
        u,cc = bools[a]
        print(f'      a{a} x_{u} coeff_in_eq={c} boolcoeff={cc} free={u in FREESET} '
              f'neqs={len(aeq[a])} nother_atoms={len(L.var_atoms[u])-1}')

print()
print('=== the -29 block (376 atoms / 347 eqs) ===')
vs = [bools[a][0] for a in S1]
print('free vars in block:', sum(1 for u in vs if u in FREESET), '/', len(vs))
print('overlap with failing eqs:', sorted(E1 & set(FAIL)))
cc = collections.Counter(bools[a][1] for a in S1)
print('bool coeff histogram in block:', cc)
h = collections.Counter(len(L.var_atoms[u])-1 for u in vs)
print('other-atoms-per-var histogram in block:', sorted(h.items()))

# decompose block into connected components (atoms sharing equations)
par = {a:a for a in S1}
def find(x):
    while par[x]!=x: par[x]=par[par[x]]; x=par[x]
    return x
def uni(x,y):
    x,y=find(x),find(y)
    if x!=y: par[x]=y
e2a = collections.defaultdict(list)
for a in S1:
    for e in aeq[a]&E1: e2a[e].append(a)
for e,As in e2a.items():
    for a in As[1:]: uni(As[0],a)
comp = collections.defaultdict(list)
for a in S1: comp[find(a)].append(a)
print(f'connected components: {len(comp)}')
rows=[]
for r,As in comp.items():
    Es = set().union(*[aeq[a] for a in As])
    rows.append((len(Es)-len(As), len(As), len(Es), sorted(As)))
rows.sort()
for defc,na,ne,As in rows[:20]:
    us=[bools[a][0] for a in As]
    print(f'  comp: {na} atoms, {ne} eqs, defic {defc}, free={sum(1 for u in us if u in FREESET)}')
    if na<=12: print(f'      atoms={As} vars={us}')

json.dump({'comps':[[r[0],r[1],r[2],r[3]] for r in rows]}, open(os.path.join(HERE,'bn_comps.json'),'w'))
