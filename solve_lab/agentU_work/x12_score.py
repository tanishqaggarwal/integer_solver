import sys, pickle, json
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentU_work')
import umodel as U, uscore as SC, checker
v0=checker.load_assignment('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json')
sd0=SC.ENG.seed_of(v0)
DRV=[642,1329,8731,9118,9413,10903,17325,18956,22162,28730,29854,30213,31864]
tgt=(v0[U.OUT[U.ROOT][0]['vab']], v0[U.OUT[U.ROOT][1]['vab']])
S={2081,24601}
rv={2081:U.LIFTC[24601], 24601:U.LIFTC[24601]}
v,isl,valn=U.assignment(S, routeval=rv, beta=U.ROOT, betaval=tgt)
s1=SC.seed_of_build(v)
print('build free-seed entries: %d'%len(s1))
n1,vv1=SC.score(s1); print('CONTROL A  build alone (no 27994 driver)      -> %d failing'%n1)
s2=SC.seed_of_build(v, {k:sd0[k] for k in DRV if k in sd0})
print('with DRV: %d entries'%len(s2))
n2,vv2=SC.score(s2); print('CONTROL B  build + deliverable DRV entries    -> %d failing'%n2)
print('   vars differing from deliverable: %d'%sum(1 for i in range(38748) if vv2[i]!=v0[i]))
# minimal DRV subset needed
only=[k for k in DRV if k not in s1]
print('DRV keys not already in build seed:',only)
