import sys, pickle, json
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentU_work/mirror')
sys.path.insert(0,'/home/user/integer_solver/solve_lab')
import harness as H
S=pickle.load(open('x_seed.pkl','rb')); sd=S['seed']; FREE=set(S['FREE'])
L=pickle.load(open('v_leaves.pkl','rb')); XY=pickle.load(open('w_xy.pkl','rb'))
D=pickle.load(open('v_defs.pkl','rb'))
sel2exp=L['sel2exp']
# reverse maps
role={}
for e,(w,m,z,C) in XY['X'].items():
    role.setdefault(w,[]).append(('Xwire',e)); role.setdefault(z,[]).append(('Xz',e))
for e,(w,m,z,C) in XY['Y'].items():
    role.setdefault(w,[]).append(('Ywire',e)); role.setdefault(z,[]).append(('Yz',e))
for s,e in sel2exp.items(): role.setdefault(s,[]).append(('SEL',e))
# z product factors
W=pickle.load(open('w_z.pkl','rb'))
print('w_z rows[0]', W['rows'][0])
for k in sorted(sd):
    r=role.get(k,[])
    print('x_%-6d bits=%-5d atoms=%-3d role=%s'%(k,sd[k].bit_length(),len(H.occ.get(k,[])),r))
print()
print('--- selectors free? ', sum(1 for s in sel2exp if s in FREE),'/256')
print('--- X wires free? ', sum(1 for e,(w,m,z,C) in XY['X'].items() if w in FREE),'/256')
print('--- Y wires free? ', sum(1 for e,(w,m,z,C) in XY['Y'].items() if w in FREE),'/256')
print('--- X z free? ', sum(1 for e,(w,m,z,C) in XY['X'].items() if z in FREE),'/256')
print('--- Y z free? ', sum(1 for e,(w,m,z,C) in XY['Y'].items() if z in FREE),'/256')
