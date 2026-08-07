import sys, json, time, pickle
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentE_work')
import engine as E, jsolve as J
from fractions import Fraction
C=J.C
base={18956:C,4279:1,26005:1}
r0,_=J.resid(base)
print("base bad", sorted(r0))
S=set()
for a in r0: S|=J.cone_free(a)
S-= {18956,4279,26005,22162,30213}
S=sorted(S)
print("stage1 free vars", len(S))
r0,cols,nonlin=J.build(base,S)
affected=set()
for f,c in cols.items(): affected|=set(c)
print("atoms affected:", len(affected), sorted(affected))
print("nonlinear pairs:", len(nonlin), sorted(nonlin)[:20])
pickle.dump({'r0':r0,'cols':cols,'nonlin':nonlin,'S':S},open('jac1.pkl','wb'))
