import os,sys; os.chdir('/home/user/integer_solver/solve_lab'); sys.path.insert(0,'.')
import heal_harness as H, json
import sy_build as B
p=H.p
# check handle vars
vA=H.loadd('best_agentA_39022.json')
for v in H.freeinp: H.val[v]=vA.get(v,0)
H.val[4287]=1; H.forward()
for name in [7927,24548,25442,5015,34661,16900]:
    print('x_%d: %s eqs, %s'%(name,'FREE' if name in H.freeinp else 'gate','val0' if H.val[name]==0 else 'nz'))
# trace x_7927
def isfree(x): return x in H.freeinp
print('x_7927 free?',isfree(7927),' x_5015 free?',isfree(5015))
