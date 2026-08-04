#!/usr/bin/env python3
import re,json
p=2**256-2**32-977; NVARS=38748
VAR_RE=re.compile(r'x_(\d+)')
FAILS=[2554,6816,8124,8680,9421,12231,12270,12350,14584,22044,29125]
codes={}
with open('../EQUATIONS.txt') as f:
    for i,line in enumerate(f):
        line=line.strip()
        if i in FAILS:
            lhs=line.rsplit('=',1)[0]
            codes[i]=compile(VAR_RE.sub(r'v[\1]',lhs),'<eq>','eval')
def loadv(path):
    d=json.load(open(path)); v=[0]*NVARS
    for k,val in d.items():
        idx=int(k[2:]) if k.startswith('x_') else int(k); v[idx]=int(val)
    return v
vA=loadv('best_agentA_39022.json')
v013=loadv('best/new_instance_partial_39013.json')
d4432=(vA[4432]-v013[4432])
d7068=(vA[7068]-v013[7068])
print(f"delta x_4432 = {d4432}  (mod p = {d4432%p})")
print(f"delta x_7068 = {d7068}  (mod p = {d7068%p})")
# coefficient of x_j in eq i = resid(v with x_j+1) - resid(v with x_j) ... but eqs are linear in these, so:
# alpha_i = d(resid)/d(x_4432): evaluate at vA and vA with x_4432 bumped by 1
def coef(i, var, base):
    ns={'v':base,'__builtins__':{}}
    r0=eval(codes[i],ns)
    base[var]+=1
    r1=eval(codes[i],ns)
    base[var]-=1
    return (r1-r0)  # linear coefficient (exact)
print("\neq :  alpha(x4432)   beta(x7068)   resid%p   check(alpha*d4432+beta*d7068)%p")
import copy
for i in FAILS:
    a=coef(i,4432,vA)%p
    b=coef(i,7068,vA)%p
    ns={'v':vA,'__builtins__':{}}
    r=eval(codes[i],ns)%p
    chk=(a*d4432+b*d7068)%p
    print(f"{i}: a={a}  b={b}  resid={r}  chk={chk}  match={r==chk}")
