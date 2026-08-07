import sys, json, collections
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentO_work')
import simO, engine as E, harness as H
d=json.load(open('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json'))
vd=[0]*E.NV
for k,x in d.items(): vd[int(k.split('_')[1])]=int(x)
bad=E.badatoms(vd)
print('bad atoms',sorted(bad))
# equation membership
eqa=getattr(H,'eqatoms',None)
print('harness attrs with eq:',[a for a in dir(H) if 'eq' in a.lower()])
