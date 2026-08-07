import json, pickle, sys, collections
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentF_work')
from fwd import Engine,NV
E=Engine()
p=115792089237316195423570985008687907853269984665640564039457584007908834671663
M=pickle.load(open('full_model.pkl','rb')); NODE=M['NODE']; OUT=M['OUT']; ROOT=M['ROOT']
D=json.load(open('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json'))
vv=[0]*NV
for k,v in D.items():
    if k.startswith('x_'): vv[int(k[2:])]=int(v)
r=E.run(vv)
nz=[i for i,x in enumerate(r) if x]
print('deliverable nonzero atoms: %d ; nonzero mod p: %d'%(len(nz),sum(1 for i in nz if r[i]%p)))
for i in nz: print('   ',E.res[i][:150])
bad=E.score(r); print('F-model failing equations:',len(bad))
# which node do these atoms belong to?
sel=NODE[ROOT]['sab'][0]
print('root sel_ab var = x%d'%sel)
d0,d1=OUT[ROOT]
print('root vab wires x%d x%d ; root vb slot wires x%d x%d'%(d0['vab'],d1['vab'],d0['vb'],d1['vb']))
