"""Minimum-weight coset decoding, step 1: the linear equations' direction classes.

The 6,614 linear equations are functionals r_e on the 6,122 unknowns.  Moving off the
common solution variety by d disturbs exactly the equations with r_e.d != 0, so the cost
in equations of ANY departure is the weight of the codeword (r_e.d)_e.  Equations whose
functionals are proportional always fail together, so the real object is the set of
DISTINCT DIRECTIONS with multiplicities.
"""
import os, sys, pickle, collections, time
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentG_work')
import gsym2 as G
from gsym2 import L, ad, P
D=pickle.load(open('/home/user/integer_solver/solve_lab/agentG_work/coset_model.pkl','rb'))
NB=D['NB']; lin=D['lin']; non=D['non']; badlin=set(D['badlin'])
n=len(NB)
rows=[]
for i,f in lin:
    r={}
    for m,c in f.items():
        if m: r[m[0][0]]=c%P
    rows.append((i,r))
cls=collections.defaultdict(list)
for i,r in rows:
    if not r: continue
    k0=min(r); iv=pow(r[k0],-1,P)
    key=tuple(sorted((c,v*iv%P) for c,v in r.items()))
    cls[key].append(i)
mult=collections.Counter(len(v) for v in cls.values())
print('linear equations %d -> distinct directions %d'%(len(rows),len(cls)))
print('multiplicity histogram:',dict(sorted(mult.items())))
lo=[(len(v),v) for v in cls.values() if len(v)<=8]
print('direction classes with multiplicity <= 8: %d'%len(lo))
for m,v in sorted(lo)[:25]: print('   mult %d : eqs %s'%(m,v[:8]))
# where do the deliverable's seven violated equations sit?
print('\ndeliverable violated equations and their classes:')
for key,v in cls.items():
    if set(v)&badlin: print('   mult %d  eqs %s   (violated: %s)'%(len(v),v[:10],sorted(set(v)&badlin)))
pickle.dump({'cls':{k:v for k,v in cls.items()},'rows':rows},open('/home/user/integer_solver/solve_lab/agentG_work/classes.pkl','wb'))
