import sys; sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentO_work')
import sweep, json
r,msg,t = sweep.run(dict(sweep.C.base))
print('baseline empty-set:',msg,'t=%.1f'%t)
if r: print('  fails',r[0],'score',39033-r[0],'bad',r[2])
