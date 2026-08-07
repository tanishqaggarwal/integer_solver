import sys, pickle, json, time
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentU_work/mirror')
sys.path.insert(0,'/home/user/integer_solver/solve_lab')
import engine3 as E3, harness as H, checker
eng=E3.Eng(E3.BASE_DEMOTE)
codes,_=checker.load_equations()
v0=checker.load_assignment('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json')
sd=eng.seed_of(v0)
base=eng.forward(sd)
f0=checker.evaluate_all(codes,base)
print('base failing',len(f0), sorted(f0)[:10])
XY=pickle.load(open('w_xy.pkl','rb'))
print('Xwire(72)=x_%d  Xconst=%d... m=%d'%(XY['X'][72][0],XY['X'][72][3]%10**20,XY['X'][72][1]))
print('seed val x_22152 == Xconst(72)? ', sd[22152]==XY['X'][72][3])
print('seed val x_33462 == Yconst(72)? ', sd[33462]==XY['Y'][72][3])
print('seed val x_6418  == Xconst(235)?', sd[6418]==XY['X'][235][3])
print('seed val x_12553 == Yconst(235)?', sd[12553]==XY['Y'][235][3])
# PERTURB the leaf-72 X wire only
for delta in (1,):
    s2=dict(sd); s2[22152]=sd[22152]+delta
    t=time.time(); v2=eng.forward(s2); el=time.time()-t
    ch=[i for i in range(38748) if v2[i]!=base[i]]
    f2=checker.evaluate_all(codes,v2)
    print('perturb x_22152 +%d: %d vars changed, %d failing (%.2fs)'%(delta,len(ch),len(f2),el))
    print('  changed vars:',ch[:40])
