import sys, collections; sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentA_work')
import env, lib as L
from amk_model import build, knobpoly, v0
P=env.P
A,K,R,rows,QUAD=build([37887,41906]),None,None,None,None
A,K,R=build([37887,41906])
kp=knobpoly(37887,K,v0)
print('knobs:',K)
print('a37887 as polynomial in the knobs (%d terms):'%len(kp))
for m,c in sorted(kp.items(), key=lambda kv:(len(kv[0]),kv[0])):
    nm='*'.join('x%d'%K[i] for i in m) or '1'
    cs=str(c) if abs(c)<10**15 else '%s(%dd)'%(str(c)[:14],len(str(abs(c))))
    print('   %-28s %s'%(nm,cs))
d0=[v0[u] for u in K]
val=0
for m,c in kp.items():
    t=c
    for i in m: t*=d0[i]
    val+=t
print('value at current knobs =',val)
