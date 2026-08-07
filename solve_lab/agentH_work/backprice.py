"""BACKWARDS step 5: price every atom of the balance-2 region that the realizable lattice
   cannot reach.  If any costs less than it buys, the backwards route beats 7."""
import frameB as FB, ev, json, time
from frameB import Frame, State
from collections import defaultdict
p=115792089237316195423570985008687907853269984665640564039457584007908834671663
D=json.load(open('backreal.json')); R=set(D['R']); S=D['S']; KN=set(int(k) for k in D['knobs'])
REACH={22229,22230,35758,35759,35760,35761,35762}
UNREACH=[a for a in S if a not in REACH]
print('atoms in the balance-2 region the lattice cannot reach:',UNREACH)
fr=Frame([642,28730,29854,31864])
W=json.load(open('../best/new_instance_partial_39026.json'))
v=[0]*38748
for k,val in W.items(): v[int(k[2:]) if k.startswith('x_') else int(k)]=int(val)
base=State(fr,{u:v[u] for u in fr.free if v[u]!=0})
B=len(base.fails); print('base failing',B)
rows=[]
for a in UNREACH:
    cands=fr.SUPV.get(a,[])
    best=None
    for X in cands:
        for d in (1,-1,p,-p):
            g=base.clone().set_free({X:base.fv.get(X,0)+d})
            if g.av[a]==base.av[a]: continue
            out=[q for q in g.av if q not in set(S) and g.av[q]!=base.av[q]]
            k=(len(g.fails),len(out))
            if best is None or k<best[0]: best=(k,X,d,len(out))
    if best is None:
        rows.append((a,None,None,None,None)); print('  a%-6d : NO move found in %d support vars'%(a,len(cands)))
    else:
        (f,o),X,d,no=best
        rows.append((a,X,d,f,no))
        print('  a%-6d : best move x_%-6d  -> failing %3d (base %d), atoms broken outside region %d'%(a,X,f,B,no))
paying=[r for r in rows if r[3] is not None and r[3]<B]
print()
print('moves that reduce the failing count below %d: %s'%(B,paying))
json.dump([[r[0],r[1],str(r[2]),r[3],r[4]] for r in rows],open('backprice.json','w'))
