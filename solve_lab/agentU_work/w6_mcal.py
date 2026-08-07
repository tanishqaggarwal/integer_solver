"""W6: CALIBRATION ONLY of M's engine as a forward propagator, before any pricing."""
import sys, json, time
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentU_work/mirror'); 
sys.path.insert(0,'/home/user/integer_solver/solve_lab')
t0=time.time()
import engine3 as E3
print('engine3 imported %.1fs ; NV=%d'%(time.time()-t0,E3.NV))
eng=E3.Eng(E3.BASE_DEMOTE)
import checker
codes,_=checker.load_equations()
v0=checker.load_assignment('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json')
print('deliverable via checker: %d failing'%len(checker.evaluate_all(codes,v0)))
sd=eng.seed_of(v0)
print('seed extracted: %d entries'%len(sd))
v1=eng.forward(sd)
n=len(checker.evaluate_all(codes,list(v1)[:38748] if len(v1)>38748 else v1))
diff=sum(1 for i in range(38748) if v0[i]!=v1[i])
print('CALIBRATION: forward(seed_of(deliverable)) -> %d failing via checker ; vars differing from deliverable: %d'%(n,diff))
print('PASS' if n==7 and diff==0 else 'FAIL')
