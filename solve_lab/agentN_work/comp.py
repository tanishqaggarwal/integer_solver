"""Frame-B region: compensator atoms and whether my extra free inputs can SET them."""
import frameB as FB, json, time, itertools
from frameB import Frame, State
p=115792089237316195423570985008687907853269984665640564039457584007908834671663
fr=Frame([642,28730,29854,31864])
W=json.load(open('../best/new_instance_partial_39026.json'))
v=[0]*38748
for k,val in W.items(): v[int(k[2:]) if k.startswith('x_') else int(k)]=int(val)
base=State(fr,{u:v[u] for u in fr.free if v[u]!=0})
NZ=set(base.nz())
E=set()
for a in NZ: E.update(fr.eq_of[a])
print('nonzero atoms',sorted(NZ)); print('region |E|=%d  satisfied %d  failing %d'%(len(E),len(E)-len(base.fails),len(base.fails)))
comp=[a for a in fr.checks if a not in NZ and fr.eq_of[a] and set(fr.eq_of[a])<=E]
print('compensator atoms (footprint entirely inside E):',len(comp),comp)
# settability: a free input that moves the atom and breaks NOTHING outside NZ|comp
allowed=NZ|set(comp)
sett={}
t0=time.time()
cands=set()
for a in list(NZ)+comp: cands.update(fr.SUPV[a])
print('free inputs reaching the region:',len(cands))
for X in cands:
    g=base.clone().set_free({X:base.fv.get(X,0)+1})
    moved=[a for a in allowed if g.av[a]!=base.av[a]]
    outside=[a for a in g.av if a not in allowed and g.av[a]!=base.av[a]]
    if moved and not outside:
        for a in moved: sett.setdefault(a,[]).append(X)
print('free inputs with ZERO collateral outside the region: atoms movable =',{a:len(x) for a,x in sett.items()})
print('  detail:',{a:x[:6] for a,x in sett.items()})
print('%.0fs'%(time.time()-t0))
# exhaustive pair scan restricted to those zero-collateral knobs
knobs=sorted({X for xs in sett.values() for X in xs})
print('zero-collateral knobs:',len(knobs),knobs[:40])
best=len(base.fails); found=None
D=[1,-1,p,-p]
t0=time.time(); n=0
for i,X in enumerate(knobs):
    for dx in D:
        gx=base.clone().set_free({X:base.fv.get(X,0)+dx})
        for Y in knobs[i+1:]:
            for dy in D:
                g=gx.clone().set_free({Y:base.fv.get(Y,0)+dy}); n+=1
                if len(g.fails)<best:
                    best=len(g.fails); found=(X,dx,Y,dy,g.clone())
                    print('  PAIR IMPROVE failing=%d  x_%d,x_%d'%(best,X,Y),flush=True)
print('pair scan over zero-collateral knobs: %d pairs, %.0fs, best failing %d'%(n,time.time()-t0,best))
if found and best<len(base.fails):
    st=found[4]; out='H_%d_pair.json'%(39033-best)
    json.dump({('x_%d'%i):st.v[i] for i in range(38748) if st.v[i]!=0},open(out,'w'))
    print('WROTE',out)
