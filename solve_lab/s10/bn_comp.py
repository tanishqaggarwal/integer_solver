"""bn_comp: components of the 311-atom maximal support; cone LP on the key-bit component."""
import os, sys, json, collections
HERE='/home/user/integer_solver/solve_lab/s10'
sys.path.insert(0, HERE)
import bn_lib as B
import lib as L, tools as T, ad
import bn_cone2 as C

bools=B.bools_map()
peel=set(json.load(open(os.path.join(HERE,'bn_cone.json')))['blk29_peel'])
inf=json.load(open(os.path.join(HERE,'bn_infl.json')))
KEYA={a for a in peel if bools[a][0] in set(inf['bool_anc2'])}
print('key atoms:',sorted(KEYA))

par={a:a for a in peel}
def find(x):
    while par[x]!=x: par[x]=par[par[x]]; x=par[x]
    return x
def uni(x,y):
    x,y=find(x),find(y)
    if x!=y: par[x]=y
e2a=collections.defaultdict(list)
for a in peel:
    for e in L.atom2eq[a]: e2a[e].append(a)
for e,As in e2a.items():
    if len(As)<2: continue
    for b in As[1:]: uni(As[0],b)
comp=collections.defaultdict(list)
for a in peel: comp[find(a)].append(a)
print('components:',len(comp))
rows=[]
for r,As in comp.items():
    E=set()
    for a in As: E|=set(L.atom2eq[a])
    fr=sum(1 for a in As if bools[a][0] in B.FREESET)
    rows.append((len(As),len(E),len(E)-len(As),fr,bool(set(As)&KEYA),sorted(As)))
rows.sort(key=lambda r:-r[0])
for na,ne,df,fr,haskey,As in rows[:12]:
    print(f'  comp {na} atoms, {ne} eqs, defic {df}, free {fr}, contains_key_bits={haskey}')
kc=[r for r in rows if r[4]]
print(f'key-bit component: {kc[0][0]} atoms, {kc[0][1]} eqs, defic {kc[0][2]}, free {kc[0][3]}')
S=kc[0][5]
json.dump({'keycomp':S}, open(os.path.join(HERE,'bn_keycomp.json'),'w'))
st,res=C.cone(set(S), tlimit=1200)
print('KEY-COMPONENT CONE RESULT:',st,flush=True)
out={'status':st,'S':S,
     'sol':({str(a):[v.numerator,v.denominator] for a,v in res.items()} if isinstance(res,dict) else None)}
json.dump(out, open(os.path.join(HERE,'bn_keycone.json'),'w'))
if isinstance(res,dict):
    for a,v in sorted(res.items()): print(f'   a{a} x_{bools[a][0]} c={bools[a][1]} t={v}')
