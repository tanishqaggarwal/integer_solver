"""BACKWARDS step 3: for the grown region, measure the REALIZABLE atom-value lattice in frame B
   and maximise the number of satisfied equations exactly over the integers."""
import frameB as FB, ev, json, time, itertools
from frameB import Frame, State
from fractions import Fraction
from collections import defaultdict
p=115792089237316195423570985008687907853269984665640564039457584007908834671663
G=json.load(open('backgrow.json'))
R=sorted(G['witness_region']); print('region |R| =',len(R))
eq_atoms=[]; atom_eqs=defaultdict(set)
for i,(m,sq,tl) in enumerate(ev.eq_terms):
    s=set(a for c,a in tl); eq_atoms.append(s)
    for a in s: atom_eqs[a].add(i)
Rset=set(R)
S=sorted(a for a in set().union(*[eq_atoms[e] for e in R]) if atom_eqs[a]<=Rset)
print('|S| =',len(S),'balance',len(R)-len(S)); print('S =',S)
fr=Frame([642,28730,29854,31864])
W=json.load(open('../best/new_instance_partial_39026.json'))
v=[0]*38748
for k,val in W.items(): v[int(k[2:]) if k.startswith('x_') else int(k)]=int(val)
base=State(fr,{u:v[u] for u in fr.free if v[u]!=0})
print('base score',base.score(),'failing in R:',len([e for e in R if e in base.fails]))
Sset=set(S)
# free inputs that move ONLY atoms of S
cands=set()
for a in S: cands.update(fr.SUPV.get(a,[]))
print('free inputs reaching S:',len(cands))
knobs={}
for X in sorted(cands):
    g=base.clone().set_free({X:base.fv.get(X,0)+1})
    diff={a:g.av[a]-base.av[a] for a in g.av if g.av[a]!=base.av[a]}
    if diff and all(a in Sset for a in diff):
        # linearity check
        h=base.clone().set_free({X:base.fv.get(X,0)+2})
        d2={a:h.av[a]-base.av[a] for a in h.av if h.av[a]!=base.av[a]}
        lin=all(d2.get(a,0)==2*diff[a] for a in diff) and set(d2)==set(diff)
        knobs[X]=(diff,lin)
print('ZERO-COLLATERAL knobs into S:',len(knobs))
for X,(d,lin) in list(knobs.items())[:40]:
    print('   x_%-6d linear=%s moves %s'%(X,lin,{a:(str(t)[:12]+'..' if len(str(t))>12 else t) for a,t in d.items()}))
json.dump({'R':R,'S':S,'knobs':{str(k):[{str(a):str(t) for a,t in d.items()},l] for k,(d,l) in knobs.items()}},
          open('backreal.json','w'))
