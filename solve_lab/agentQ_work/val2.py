import sys, os, json
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
E.run(v)
Lv=leaves(); pins=json.load(open(os.path.join(HERE,'..','agentF_work','pins.json')))
on=[int(g) for g in pins if v[int(g)]%p==1]
tot=None
for g in on: tot=add(tot,Lv[g])
rx,ry=(tot[0]-cs)%p, tot[1]%p
hit={}
for i in range(NV):
    m=v[i]%p
    if m==rx: hit.setdefault('sum_x',[]).append(i)
    if m==ry: hit.setdefault('sum_y',[]).append(i)
    if m==(Lv[on[0]][0]-cs)%p: hit.setdefault('leaf%d_x'%on[0],[]).append(i)
    if m==(Lv[on[1]][0]-cs)%p: hit.setdefault('leaf%d_x'%on[1],[]).append(i)
for k,vv in hit.items(): print(k, len(vv), vv[:12])
# nonzero atoms
nz=[a for a in E.res if E.evalatom(a,v)%1==0] if hasattr(E,'evalatom') else None
print('engine attrs:', [a for a in dir(E) if not a.startswith('_')][:30])
