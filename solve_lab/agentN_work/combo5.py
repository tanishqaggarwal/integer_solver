"""Rank-raising sweep, stage 3: search integer combinations of the five a22231-movers for one
   that leaves EVERY atom outside the 8-atom region at zero.  Ground truth by direct evaluation."""
import frameB as FB, ev, json, time, itertools
from frameB import Frame, State
p=115792089237316195423570985008687907853269984665640564039457584007908834671663
fr=Frame([642,28730,29854,31864])
W=json.load(open('../best/new_instance_partial_39026.json'))
v=[0]*38748
for k,val in W.items(): v[int(k[2:]) if k.startswith('x_') else int(k)]=int(val)
base=State(fr,{u:v[u] for u in fr.free if v[u]!=0})
IN8={22229,22230,22231,35758,35759,35760,35761,35762}
MOV=[2081,4287,4432,12553,28730]
print('base failing',len(base.fails))
# individual signatures
for X in MOV:
    for d in (1,-1,p):
        g=base.clone().set_free({X:base.fv.get(X,0)+d})
        out=sorted(a for a in g.av if g.av[a]!=base.av[a] and a not in IN8)
        print('  x_%-6d d=%-4s outside-region atoms moved: %d  %s  failing=%d'%(
            X,'p' if d==p else d,len(out),out[:8],len(g.fails)))
t0=time.time(); best=(len(base.fails),None); tested=0
RNG=[-3,-2,-1,0,1,2,3]
found=[]
for c in itertools.product(RNG,repeat=5):
    if all(x==0 for x in c): continue
    tested+=1
    g=base.clone()
    ch={}
    for X,k in zip(MOV,c):
        if k: ch[X]=base.fv.get(X,0)+k
    g.set_free(ch)
    if 22231 not in [a for a in g.av if g.av[a]!=base.av[a]]: continue
    out=[a for a in g.av if g.av[a]!=base.av[a] and a not in IN8]
    if not out:
        found.append(c); print('  ZERO-COLLATERAL COMBO',c,'failing',len(g.fails),flush=True)
    if len(g.fails)<best[0]:
        best=(len(g.fails),g.clone()); print('  IMPROVE failing=%d combo=%s'%(best[0],c),flush=True)
print('combinations tested %d in %.0fs; zero-collateral combos: %d; best failing %d'%(tested,time.time()-t0,len(found),best[0]))
# also try p-scaled combinations on the two cheapest movers
t0=time.time()
for a in range(-3,4):
    for b in range(-3,4):
        for s in (1,p):
            for t in (1,p):
                if a==0 and b==0: continue
                g=base.clone().set_free({28730:base.fv.get(28730,0)+a*s,12553:base.fv.get(12553,0)+b*t})
                out=[q for q in g.av if g.av[q]!=base.av[q] and q not in IN8]
                if not out: print('  ZERO-COLLATERAL (28730,12553)=(%d*%s,%d*%s) failing %d'%(a,'p' if s==p else 1,b,'p' if t==p else 1,len(g.fails)),flush=True)
                if len(g.fails)<best[0]:
                    best=(len(g.fails),g.clone()); print('  IMPROVE failing=%d'%best[0],flush=True)
print('p-scaled pair sweep %.0fs, best failing %d'%(time.time()-t0,best[0]))
if best[0]<len(base.fails):
    st=best[1]; sc=39033-best[0]
    json.dump({('x_%d'%i):st.v[i] for i in range(38748) if st.v[i]!=0},open('H_%d_rank.json'%sc,'w'))
    print('WROTE H_%d_rank.json'%sc)
else:
    print('NO combination of the five a22231-movers beats 7.')
