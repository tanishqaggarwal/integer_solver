import sys, pickle
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentU_work')
sys.path.insert(0,'/home/user/integer_solver/solve_lab')
import umodel as U, checker
v0=checker.load_assignment('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json')
print('SLOTS (both sides live):',len(U.SLOTS),' ROOT in SLOTS:',U.ROOT in U.SLOTS)
S={2081,24601}
# deliverable: route carries leaf-72's lifted point on BOTH sides; pins honest; beta=ROOT driven to target
TGT=(v0[U.OUT[U.ROOT][0]['vab']], v0[U.OUT[U.ROOT][1]['vab']])
print('TARGET =',[str(x)[:40] for x in TGT], [x.bit_length() for x in TGT])
rv={2081:U.LIFT[24601], 24601:U.LIFT[24601]}
v,isl,valn=U.assignment(S, routeval=rv, beta=U.ROOT, betaval=TGT)
print('built %d vars'%len(v))
# compare against deliverable on the vars we set
diff=[k for k in v if v0[k]!=v[k]]
print('vars where my build differs from deliverable: %d'%len(diff))
for k in diff[:20]: print('   x_%d mine=%s deliv=%s'%(k,str(v[k])[:30],str(v0[k])[:30]))
# deliverable nonzero vars not in my build
extra=[k for k in range(38748) if v0[k]!=0 and k not in v]
print('deliverable nonzero vars NOT in my build: %d'%len(extra))
print('   ',extra[:25])
