import sys, time, json, os; sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentO_work')
import simO
for bit in [22492,13710]:
    s=dict(simO.C.base); s[bit]=1
    t=time.time()
    n,ns,av,v,blk,info=simO.solve_once(s,{bit},maxr=5,maxv=4000,verbose=True)
    print(bit,'fails',n,'score',39033-n,'bad',av,'info',info,'t=%.0f'%(time.time()-t),flush=True)
