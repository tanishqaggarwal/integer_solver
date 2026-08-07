import sys, json; sys.path.insert(0,'.')
from fractions import Fraction as F
from math import gcd
M=json.load(open('/home/user/integer_solver/solve_lab/agentA_work/model22.json'))
K=M['K']; rows=M['rows']
aff=[r for r in rows if not r['quad']]
Lin=[[F(int(x)) for x in r['lin']] for r in aff]
C=[F(int(r['c'])) for r in aff]
EQ=[r['eq'] for r in aff]
n=len(aff); m=len(K)
# solve Lin*d = -C exactly (overdetermined, consistent, rank m)
aug=[Lin[i]+[-C[i]] for i in range(n)]
piv=[];r=0
for c in range(m):
    pr=None
    for i in range(r,n):
        if aug[i][c]!=0: pr=i;break
    if pr is None: continue
    aug[r],aug[pr]=aug[pr],aug[r]
    pv=aug[r][c]; aug[r]=[x/pv for x in aug[r]]
    for i in range(n):
        if i!=r and aug[i][c]!=0:
            f=aug[i][c]; aug[i]=[a-f*b for a,b in zip(aug[i],aug[r])]
    piv.append(c); r+=1
sol=[F(0)]*m
for i,c in enumerate(piv): sol[c]=aug[i][m]
# check consistency
bad=[EQ[i] for i in range(n) if sum(Lin[i][j]*sol[j] for j in range(m))+C[i]!=0]
print('consistent over Q:', not bad, bad)
d0=[int(x) for x in M['d0']]
print('%-8s %-6s %s'%('knob','denom','solution (or delta from current)'))
allint=True
for j,u in enumerate(K):
    den=sol[j].denominator
    if den!=1: allint=False
    dd=sol[j]-d0[j]
    s=str(sol[j].numerator)
    print('x%-7d den=%-30s val=%s%s  delta=%s'%(u,den if den<10**12 else '%dd:%s'%(len(str(den)),str(den)[:20]),
          s[:30]+('...' if len(s)>30 else ''),'' ,str(dd.numerator)[:24]+'/'+str(dd.denominator)[:24]))
print('ALL INTEGRAL:',allint)
json.dump({'K':K,'sol_num':[str(x.numerator) for x in sol],'sol_den':[str(x.denominator) for x in sol]},
          open('/home/user/integer_solver/solve_lab/agentA_work/qsol.json','w'))
