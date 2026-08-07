"""Rank-raising sweep, stage 1: collateral census over ALL free inputs of frame B.
   CROSSOVER ARGUMENT: adding one knob raises the lattice rank by at most 1, so the number of
   zeroable rows rises by at most 1.  Every NEW equation dragged in adds at least 1 to |R'|.
   So a knob can only pay if it drags in ZERO new equations, i.e. every atom it moves has its
   whole equation footprint inside the 12-equation region.  This makes the sweep exhaustive."""
import frameB as FB, ev, json, time
from frameB import Frame, State
from collections import defaultdict
p=115792089237316195423570985008687907853269984665640564039457584007908834671663
atom_eqs=defaultdict(set)
for i,(m,sq,tl) in enumerate(ev.eq_terms):
    for c,a in tl: atom_eqs[a].add(i)
fr=Frame([642,28730,29854,31864])
W=json.load(open('../best/new_instance_partial_39026.json'))
v=[0]*38748
for k,val in W.items(): v[int(k[2:]) if k.startswith('x_') else int(k)]=int(val)
base=State(fr,{u:v[u] for u in fr.free if v[u]!=0})
BASE7={22229,22230,35758,35759,35760,35761,35762}
R=set()
for a in BASE7: R|=atom_eqs[a]
print('base region |R| =',len(R),'failing',len(base.fails))
INSIDE=set(a for a in fr.checks if atom_eqs[a] and atom_eqs[a]<=R)
print('atoms with footprint entirely inside R:',sorted(INSIDE))
t0=time.time()
census=defaultdict(list)          # newly-dragged-in equation count -> knobs
moves22231=[]
DELTAS=[1,p]
for i,X in enumerate(fr.free):
    cur=base.fv.get(X,0)
    moved=set()
    for d in DELTAS:
        g=base.clone().set_free({X:cur+d})
        moved |= set(a for a in g.av if g.av[a]!=base.av[a])
    if not moved: continue
    newR=set()
    for a in moved: newR|=atom_eqs[a]
    extra=len(newR-R)
    census[extra].append(X)
    if 22231 in moved: moves22231.append((X,extra,sorted(moved)))
    if i%2000==0: print('  %d/%d %.0fs'%(i,len(fr.free),time.time()-t0),flush=True)
print('\ncollateral census (new equations dragged in -> #free inputs):')
for k in sorted(census)[:12]:
    print('   +%d new eqs : %d knobs  %s'%(k,len(census[k]),census[k][:12]))
print('\nknobs that move a22231 (the only compensator inside R): %d'%len(moves22231))
for X,e,mv in moves22231[:20]:
    print('   x_%-6d drags in %2d new equations, moves atoms %s'%(X,e,mv))
json.dump({'census':{str(k):v for k,v in census.items()},
           'moves22231':[[X,e,mv] for X,e,mv in moves22231]},open('knobcensus.json','w'))
print('\nCROSSOVER: only knobs with +0 new equations can possibly pay.  Those are: %s'%census.get(0,[]))
