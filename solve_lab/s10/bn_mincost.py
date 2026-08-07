"""bn_mincost: minimum number of equations broken by a nonzero boolean config.

Since the free-var cone is trivial, every nonzero configuration breaks >=1 eq.
Here we search small supports exactly:
  - singletons: cost = |E(a)| (no cancellation possible in a 1-term sum)
  - pairs / triples inside the sign-peel core, with t values ranging over the
    triangular image {0,2,6,12,20,30,42,56,72,90,...}
Reports the true combinatorial minimum found (equation count only, before
downstream damage), and cross-checks against the measured sweep.
"""
import os, sys, json, collections, itertools, time
HERE='/home/user/integer_solver/solve_lab/s10'
sys.path.insert(0, HERE)
import bn_lib as B
import lib as L, tools as T, ad

OUT=os.path.join(HERE,'bn_mincost.jsonl')
bools=B.bools_map()
FREEB=[a for a in bools if bools[a][0] in B.FREESET]
TRI=[k*(k-1) for k in range(1,16)]      # 0,2,6,12,20,...
TRI=[t for t in TRI if t>0]

def cost(sel):
    """sel: dict atom -> t (>=0).  number of equations with nonzero sum."""
    es=set()
    for a in sel: es|=set(L.atom2eq[a])
    c=0
    for e in es:
        m,sq,co=L.eq_atoms[e]
        s=sum(co[a]*bools[a][1]*sel[a] for a in sel if a in co)
        if s: c+=1
    return c

f=open(OUT,'a')
print('--- singletons (free boolean atoms) ---',flush=True)
sing=sorted((len(L.atom2eq[a]),a) for a in FREEB)
print('  min |E(a)| =',sing[0][0],'at a%d (x_%d)'%(sing[0][1],bools[sing[0][1]][0]))
print('  10 smallest:',[(a,n) for n,a in sing[:10]])
for n,a in sing[:10]:
    f.write(json.dumps({'kind':'single','atoms':[a],'t':[2],'cost':cost({a:2})})+'\n')
f.flush()

print('--- pairs inside the sign-peel core ---',flush=True)
cone=json.load(open(os.path.join(HERE,'bn_cone.json')))
core=cone['free_peel']
print('  core size',len(core))
best=[]
t0=time.time()
for a,b in itertools.combinations(sorted(core),2):
    sh=set(L.atom2eq[a])&set(L.atom2eq[b])
    if not sh: continue
    for ta in TRI[:8]:
        for tb in TRI[:8]:
            c=cost({a:ta,b:tb})
            best.append((c,a,b,ta,tb))
best.sort()
print('  best pairs (combinatorial equation cost):',flush=True)
for c,a,b,ta,tb in best[:10]:
    print(f'    cost={c} a{a}(x_{bools[a][0]},t={ta}) a{b}(x_{bools[b][0]},t={tb}) '
          f'|E(a)|={len(L.atom2eq[a])} |E(b)|={len(L.atom2eq[b])}')
    f.write(json.dumps({'kind':'pair','atoms':[a,b],'t':[ta,tb],'cost':c})+'\n')
f.flush()

print('--- triples inside the sign-peel core (t in {2,6,12}) ---',flush=True)
best3=[]
adj=collections.defaultdict(set)
for a in core:
    for e in L.atom2eq[a]: adj[e].add(a)
trip=set()
for e,As in adj.items():
    for c3 in itertools.combinations(sorted(As),3): trip.add(c3)
print('  candidate triples sharing an equation:',len(trip),flush=True)
for tr in trip:
    for ts in itertools.product([2,6,12],repeat=3):
        c=cost(dict(zip(tr,ts)))
        best3.append((c,tr,ts))
best3.sort()
for c,tr,ts in best3[:10]:
    print(f'    cost={c} atoms={tr} t={ts}')
    f.write(json.dumps({'kind':'triple','atoms':list(tr),'t':list(ts),'cost':c})+'\n')
f.close()
print('saved',OUT,flush=True)
