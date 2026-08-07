import sys, pickle
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentU_work')
sys.path.insert(0,'/home/user/integer_solver/solve_lab')
import umodel as U, checker
v0=checker.load_assignment('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json')
S={2081,24601}
TGT={'X':None,'Y':None}
# ROOT axis order
print('AXIS[ROOT-ish]: OUT[ROOT] vab wires',[d['vab'] for d in U.OUT[U.ROOT]])
tgt=(v0[U.OUT[U.ROOT][0]['vab']], v0[U.OUT[U.ROOT][1]['vab']])
rv={2081:U.LIFTC[24601], 24601:U.LIFTC[24601]}
v,isl,valn=U.assignment(S, routeval=rv, beta=U.ROOT, betaval=tgt)
diff=[k for k in v if v0[k]!=v[k]]
print('built %d vars ; differing from deliverable: %d'%(len(v),len(diff)))
for k in diff: print('   x_%d mine=%s deliv=%s'%(k,str(v[k])[:34],str(v0[k])[:34]))
pickle.dump({'tgt':tgt},open('x_tgt.pkl','wb'))
