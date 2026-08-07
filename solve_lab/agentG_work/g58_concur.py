"""Given a 2-variable departure, find the point minimising the number of nonzero
equations: intersect every pair of the varying curves and count how many vanish there."""
import os, sys, pickle, itertools, collections
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentG_work')
import gsym2 as G
from gsym2 import L, ad, P
import flint
d=pickle.load(open(sys.argv[1],'rb')); polys=d['polys']; SUP=d['SUP']
var=[(i,g) for i,g in polys if any(any(e) for e in g)]
print('varying equations:',len(var),'support',SUP)
def ev(g,a,b):
    s=0
    for m,c in g.items():
        s=(s+c*pow(a,m[0],P)*pow(b,m[1],P))%P
    return s
def count(a,b): return [i for i,g in var if ev(g,a,b)]
print('at d=0: failing %d'%len(count(0,0)))
aff=[(i,g) for i,g in var if max(sum(m) for m in g)<=1]
print('affine curves: %d ; higher degree: %s'%(len(aff),[i for i,g in var if max(sum(m) for m in g)>1]))
cands=set([(0,0)])
for (i1,g1),(i2,g2) in itertools.combinations(aff,2):
    a1=g1.get((1,0),0); b1=g1.get((0,1),0); c1=g1.get((0,0),0)
    a2=g2.get((1,0),0); b2=g2.get((0,1),0); c2=g2.get((0,0),0)
    det=(a1*b2-a2*b1)%P
    if det==0: continue
    iv=pow(det,-1,P)
    x=((-c1)*b2-(-c2)*b1)*iv%P
    y=(a1*(-c2)-a2*(-c1))*iv%P
    cands.add((x,y))
print('candidate points:',len(cands))
best=None
tally=collections.Counter()
for (a,b) in cands:
    f=count(a,b); tally[len(f)]+=1
    if best is None or len(f)<len(best[1]): best=((a,b),f)
print('failing-count histogram over candidates:',dict(sorted(tally.items())))
print('\nBEST: %d failing equations -> exact mod-p score %d'%(len(best[1]),L.NEQ-len(best[1])))
print('  point d =',best[0])
print('  failing:',best[1])
pickle.dump({'point':best[0],'fail':best[1],'SUP':SUP},open(sys.argv[1].replace('dep_','best_'),'wb'))
