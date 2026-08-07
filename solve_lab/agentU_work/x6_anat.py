import sys, pickle, time
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentU_work/mirror')
sys.path.insert(0,'/home/user/integer_solver/solve_lab')
import engine3 as E3, harness as H, checker
eng=E3.Eng(E3.BASE_DEMOTE); codes,_=checker.load_equations()
v0=checker.load_assignment('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json')
sd=eng.seed_of(v0)
XY=pickle.load(open('w_xy.pkl','rb'))
CX=XY['X'][72][3]; CY=XY['Y'][72][3]
RX=sorted(k for k in sd if sd[k]==CX); RY=sorted(k for k in sd if sd[k]==CY)
SEL=[2081,24601]; LW=[6418,12553,22152,33462]
DRV=sorted(set(sd)-set(RX)-set(RY)-set(SEL)-set(LW))
def sc(s):
    return len(checker.evaluate_all(codes,eng.forward(s)))
print('groups: RX=%d RY=%d SEL=%d LW=%d DRV=%d'%(len(RX),len(RY),len(SEL),len(LW),len(DRV)))
print('DRV =',DRV)
print('full                 ->',sc(sd))
for name,g in [('RX',RX),('RY',RY),('DRV',DRV),('SEL',SEL),('RX+RY',RX+RY)]:
    s={k:v for k,v in sd.items() if k not in set(g)}
    print('minus %-8s        -> %d'%(name,sc(s)))
print('only SEL+LW          ->',sc({k:sd[k] for k in SEL+LW}))
print('only SEL+LW+RX+RY    ->',sc({k:sd[k] for k in SEL+LW+RX+RY}))
print('only SEL             ->',sc({k:sd[k] for k in SEL}))
print('empty                ->',sc({}))
print()
print('leave-one-out:')
for k in sorted(sd):
    s=dict(sd); del s[k]
    print('  drop x_%-6d -> %d'%(k,sc(s)))
