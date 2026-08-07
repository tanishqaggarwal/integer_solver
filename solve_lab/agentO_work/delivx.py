import sys, json
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentO_work')
import simO, engine as E, harness as H
d=json.load(open('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json'))
v=[0]*E.NV
for k,x in d.items(): v[int(k.split('_')[1])]=int(x)
bad=E.badatoms(v); ff=E.eqfails(bad)
print('deliverable: bad atoms',sorted(bad),'fails',len(ff))
print('x_7715',v[7715],'x_34554',v[34554],'x_15298',v[15298])
print('x_11150 bits',abs(v[11150]).bit_length(),'x_25739',abs(v[25739]).bit_length(),'x_37758',abs(v[37758]).bit_length())
for a in sorted(bad):
    print(' ',a,H.atoms[a][:100],'| resid bits',abs(bad[a]).bit_length())
# selector bits on?
import channels as C
print('cluster bools on:',[f for f in C.BOOLS if v[f]!=0][:40])
