import sys, os, json, collections, pickle
F='/home/user/integer_solver/solve_lab/agentV_work/mirror/F'; sys.path.insert(0,F)
from fwd import Engine, NV
from parse import node_str
from circ2 import vars_of
E=Engine()
p=115792089237316195423570985008687907853269984665640564039457584007908834671663
A=json.load(open('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json'))
v=[0]*NV
for k,x in A.items(): v[int(k[2:])]=int(x)
r=E.run(v)
bad=E.score(r)
print('failing eq idx',bad)
nz=[i for i,x in enumerate(r) if x]
print('nonzero residual atoms',len(nz))
for i in nz:
    a=E.res[i]
    print('  ',a[:250])
    print('     value mod p ==0 ?', r[i]%p==0, ' bits', r[i].bit_length())
