import json, pickle, sys
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentV_work/mirror/F')
from fwd import Engine,NV
E=Engine()
p=115792089237316195423570985008687907853269984665640564039457584007908834671663
M=pickle.load(open('full_model.pkl','rb')); NODE=M['NODE']; OUT=M['OUT']; ROOT=M['ROOT']
D=json.load(open('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json'))
g=lambda i:int(D.get('x_%d'%i,0))
print('parent slot wires  x4432=%s'%(g(4432)%p))
print('                   x7068=%s'%(g(7068)%p))
print('child vab wires    x8731=%s'%(g(8731)%p))
print('                   x9118=%s'%(g(9118)%p))
print('slot == vab (coord-wise, some order)?',
      {g(4432)%p,g(7068)%p}=={g(8731)%p,g(9118)%p})
# locate the site in my model
site=None
for n in NODE:
    for s in ('va','vb'):
        if [d[s] for d in OUT[n]]==[4432,7068] or [d[s] for d in OUT[n]]==[7068,4432]:
            site=(n,s)
print('parent node/side for wires (4432,7068):',site)
c = NODE[site[0]]['a'] if site[1]=='va' else NODE[site[0]]['b']
print('child node c = x%d ; its vab wires = %s'%(c,[d['vab'] for d in OUT[c]]))
print('c sel_ab var x%d value=%s'%(NODE[c]['sab'][0], g(NODE[c]['sab'][0])))
