import sys, json, pickle, time
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentU_work/mirror')
sys.path.insert(0,'/home/user/integer_solver/solve_lab')
import engine3 as E3, harness as H
eng=E3.Eng(E3.BASE_DEMOTE)
import checker
v0=checker.load_assignment('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json')
sd=eng.seed_of(v0)
print('FREE size', len(eng.FREE), ' SEQ', len(eng.SEQ))
print('nonzero seed entries', len(sd))
for k in sorted(sd):
    val=sd[k]
    print('x_%d  bits=%d  %s'%(k, val.bit_length(), str(val)[:60]))
pickle.dump({'FREE':eng.FREE,'seed':sd}, open('x_seed.pkl','wb'))
