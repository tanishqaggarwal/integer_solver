import os,sys; os.chdir('/home/user/integer_solver/solve_lab'); sys.path.insert(0,'.')
import heal_harness as H, re
from collections import Counter
vA=H.loadd('best_agentA_39022.json')
for v in H.freeinp: H.val[v]=vA.get(v,0)
H.val[4287]=1; H.forward()
eqcount=Counter(); VAR=re.compile(r'x_(\d+)')
for L in open('../EQUATIONS.txt'):
    if not L.strip(): continue
    for m in set(VAR.findall(L)): eqcount[int(m)]+=1
for name in [9106,9629,2239,23754,31731,35619,11368,24559,32491,27676,7574,21279]:
    tag = 'FREE' if name in H.freeinp else 'gate'
    print('x_%d: val=%s  %s  eqs=%d'%(name,H.val[name],tag,eqcount[name]))
