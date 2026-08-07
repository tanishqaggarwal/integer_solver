import sys, os, json, re
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE)
sys.path.insert(0, os.path.join(HERE,'..','agentF_work'))
from qgrp import *
from fwd import Engine, NV
E=Engine()
asg=json.load(open(os.path.join(HERE,'..','best','new_instance_partial_39026.json')))
v=[0]*NV
for k,x in asg.items():
    i=int(k.split('_')[1]) if '_' in k else int(k)
    if i<NV: v[i]=int(x)
r=E.run(v)
nz=[(E.res[i], r[i]) for i,x in enumerate(r) if x]
bad=E.score(r)
print('score', 39033-len(bad), 'failing eqs', len(bad))
print('nonzero residual atoms:', len(nz))
for a,val in nz:
    print('  ', a[:150])
    print('      value mod p =', val%p, ' | val/p integral?', val%p==0)
