import sys
sys.path.insert(0,'/home/user/integer_solver/solve_lab/s9/pins')
from build import *
codes,_=H.load_equations()
v0=list(BASE)
print('BASE residual',nz(v0),'defects',defects(v0))
print('  fails',len(H.evaluate(codes,v0)))
v=build([2081,24601],close=False)
print('close=False residual',nz(v),'defects',defects(v))
f=H.evaluate(codes,v); print('  fails',len(f),f)
# which vars differ from BASE
diff=[i for i in range(NV) if v[i]!=BASE[i]]
print('vars changed vs BASE:',len(diff))
