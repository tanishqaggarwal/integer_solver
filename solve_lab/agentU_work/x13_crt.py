import sys, pickle, math, time
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentU_work')
import umodel as U, ucrt as CR, uscore as SC, checker
v0=checker.load_assignment('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json')
tgt=(v0[U.OUT[U.ROOT][0]['vab']], v0[U.OUT[U.ROOT][1]['vab']])
sd0=SC.ENG.seed_of(v0)
DRV=[642,1329,8731,9118,9413,10903,17325,18956,22162,28730,29854,31864]
DRVSEED={k:sd0[k] for k in DRV if k in sd0}
# sanity: chain_to_root at ROOT is identity
print('chain_to_root(ROOT)',CR.chain_to_root(U.ROOT))
print('betaval_for(ROOT,tgt)==tgt:',CR.betaval_for(U.ROOT,tgt)==tgt)
# free-factor check
FREE=SC.FREE
nf=sum(1 for w,(z,m,fa,fb,ffa,ffb) in CR.ZINFO.items() if fa in FREE and fb in FREE)
print('pin wires whose BOTH z-factors are FREE in M engine: %d/%d'%(nf,len(CR.ZINFO)))
# ---- CRT construction at the ROOT slot, deliverable's own leaf pair
a,b=24601,2081
W=CR.pair_ok(a,b); print('CRT W for (72,235):', 'None' if W is None else str(W)[:40], W.bit_length() if W else '')
if W:
    Wy=U.LIFTC[a]['Y']
    ex={}
    for s in (a,b):
        e=CR.leaf_extras(s,W,Wy)
        print('  leaf_extras sel %d ->'%s, 'None' if e is None else {k:str(vv)[:18] for k,vv in e.items()})
        if e: ex.update(e)
    rv={a:{'X':W,'Y':Wy}, b:{'X':W,'Y':Wy}}
    v,isl,valn=U.assignment({a,b},routeval=rv,pinval=rv,beta=U.ROOT,betaval=tgt)
    s=SC.seed_of_build(v); s.update({k:val for k,val in ex.items() if val!=0})
    n,_=SC.score(s); print('CRT @ROOT, pair(72,235), no DRV  -> %d failing (seed %d)'%(n,len(s)))
    s2=dict(s); s2.update(DRVSEED)
    n2,_=SC.score(s2); print('CRT @ROOT, pair(72,235), + DRV  -> %d failing'%n2)
