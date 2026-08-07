"""Rank-raising sweep, stage 4 (final): a37887 = Q^2 lives in ONE equation.  x_28730 is the unique
   minimal a22231-mover and its only outside atom is a37887.  Is there ANY free input that moves
   a37887 so the two can cancel while a22231 survives?"""
import frameB as FB, ev, json, time
from frameB import Frame, State
from collections import defaultdict
p=115792089237316195423570985008687907853269984665640564039457584007908834671663
atom_eqs=defaultdict(set)
for i,(m,sq,tl) in enumerate(ev.eq_terms):
    for c,a in tl: atom_eqs[a].add(i)
print('a37887 lives in equations:',sorted(atom_eqs[37887]))
print('a22231 lives in equations:',sorted(atom_eqs[22231]))
fr=Frame([642,28730,29854,31864])
W=json.load(open('../best/new_instance_partial_39026.json'))
v=[0]*38748
for k,val in W.items(): v[int(k[2:]) if k.startswith('x_') else int(k)]=int(val)
base=State(fr,{u:v[u] for u in fr.free if v[u]!=0})
IN8={22229,22230,22231,35758,35759,35760,35761,35762}
t0=time.time(); movers=[]
for i,X in enumerate(fr.free):
    cur=base.fv.get(X,0)
    for d in (1,-1,p):
        g=base.clone().set_free({X:cur+d})
        if g.av[37887]!=base.av[37887]:
            mv=sorted(a for a in g.av if g.av[a]!=base.av[a])
            movers.append((X,d if d!=p else 'p',mv,len(g.fails)))
            break
    if i%2500==0: print('  %d/%d %.0fs'%(i,len(fr.free),time.time()-t0),flush=True)
print('\nfree inputs that move a37887: %d'%len(movers))
clean=[m for m in movers if all(a in IN8 or a==37887 for a in m[2])]
print('...of which move NOTHING outside the region except a37887: %d'%len(clean))
for X,d,mv,f in clean[:20]:
    print('   x_%-6d d=%-3s moves %s  failing=%d  (moves a22231: %s)'%(X,d,mv,f,22231 in mv))
withA=[m for m in clean if 22231 in m[2]]
without=[m for m in clean if 22231 not in m[2]]
print('\n  clean movers that DO move a22231 : %d  %s'%(len(withA),[m[0] for m in withA]))
print('  clean movers that do NOT move a22231 (the cancellers we need): %d  %s'%(len(without),[m[0] for m in without]))
if not without:
    print('\n=> a37887 cannot be cancelled.  Every direction that moves a22231 also moves a37887,')
    print('   and a37887 lives in exactly one equation, eq %s, which is therefore forced to fail.'%sorted(atom_eqs[37887]))
    print('   a22231 buys at most 1 row and costs exactly that row.  **7 IS THE FLOOR.**')
json.dump({'a37887_eqs':sorted(atom_eqs[37887]),'a22231_eqs':sorted(atom_eqs[22231]),
           'n_movers':len(movers),'clean':[[m[0],str(m[1]),m[2],m[3]] for m in clean]},open('q37887.json','w'))
