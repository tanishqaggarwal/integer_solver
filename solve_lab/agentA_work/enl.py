"""Explicit enlargement test on the 39,026 region using regsolve2 (strictly linear)."""
import sys, json, collections; sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentA_work')
from fractions import Fraction as F
import env, lib as L
from regsolve2 import build, qsolve
P=env.P
v=L.load('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json')
E=sorted(set(e for a in [22229,22230,35758,35759,35760,35761,35762] for e in L.atom2eq[a]))
A0=set(a for e in E for a in L.eq_atoms[e][2])
X9118=[1465,8263,36088,40005,40121]
X8731=[1459,8261,40005]
X7068=[2202,16897,21113,38521,39166,40066,40932]
X2099=[29090]
X4432=[2200,21114,32910,37887,39166,40066,40932]
X19964=[1461,37887,40005]
TESTS=[('base',[]),
       ('+a37887,a41906',[37887,41906]),
       ('+x9118',[37887,41906]+X9118),
       ('+x8731',[37887,41906]+X8731),
       ('+x9118+x8731',[37887,41906]+X9118+X8731),
       ('+x2099',[37887,41906]+X2099),
       ('+x7068',[37887,41906]+X7068),
       ('+x4432+x19964',[37887,41906]+X4432+X19964),
       ('+ALL',[37887,41906]+X9118+X8731+X7068+X2099+X4432+X19964),
       ]
for name,ex in TESTS:
    A=A0|set(ex)
    K,R,rows=build(v,A)
    sol,free,incons,r,piv,mat,eqs=qsolve(rows,len(K))
    skipped=sum(1 for x in rows if x[3])
    bad=[]
    if not incons:
        for j,u in enumerate(K):
            if sol[j] is None: continue
            d=sol[j].denominator
            if d!=1: bad.append((u,'p' if d==P else ('p*%d'%(d//P) if d%P==0 else str(d))))
    inK=lambda u: u in set(K)
    print('%-16s atoms=%-4d knobs=%-4d eqs=%-5d skip=%-3d rank=%-4d free=%-3d incons=%-3d nonint=%-2d %s'%(
        name,len(A),len(K),len(R),skipped,r,len(free),len(incons),len(bad),bad[:7]),flush=True)
    print('                 x9118knob=%s x8731knob=%s x2099knob=%s x7068knob=%s x4432knob=%s x19964knob=%s'%(
        inK(9118),inK(8731),inK(2099),inK(7068),inK(4432),inK(19964)))
