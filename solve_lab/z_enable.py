import os,sys,json
os.chdir('/home/user/integer_solver/solve_lab'); sys.path.insert(0,'.')
import heal_harness as H
p=H.p
vA=H.loadd('best_agentA_39022.json')
for x in H.freeinp: H.val[x]=vA.get(x,0)
H.forward()
for v in [15298,1222,35723,28505,32083,7715,34554]:
    st='FREE' if v in H.freeinp else 'gate'
    rhs=''
    if v in H.definer:
        gi=H.definer[v]; rhs=H.gates[gi][1][:50]
    print(f"x_{v}={H.val[v]%p}  {st} rhs={rhs}")
