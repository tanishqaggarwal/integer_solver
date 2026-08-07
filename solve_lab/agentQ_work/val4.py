import sys, os, json
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE)
sys.path.insert(0, os.path.join(HERE,'..','agentF_work'))
from qgrp import *
from fwd import Engine, NV, compile_node
E=Engine()
asg=json.load(open(os.path.join(HERE,'..','best','new_instance_partial_39026.json')))
v=[0]*NV
for k,x in asg.items():
    i=int(k.split('_')[1]) if '_' in k else int(k)
    if i<NV: v[i]=int(x)
# evaluate ALL atoms (definition + residual) at the assignment AS GIVEN
allatoms=sorted(set(E.atoms))
src='r[:]=['+','.join(compile_node(E.atoms[a]) for a in allatoms)+']'
prog=compile(src,'<all>','exec'); r=[0]*len(allatoms)
exec(prog,{'v':v,'r':r,'__builtins__':{}})
idx={a:i for i,a in enumerate(allatoms)}
nz=[(a,r[idx[a]]) for a in allatoms if r[idx[a]]]
print('nonzero atoms at the assignment as given:', len(nz))
for a,val in nz: print('   %-90s  |val| bits %d  mod p = %s' % (a[:90], val.bit_length(), val%p))
bad=[]
for i,row in enumerate(E.eqrows):
    t=sum(k*r[idx[a]] for k,a in row)
    if t: bad.append(i)
print('failing equation rows:', len(bad), '=> score', 39033-len(bad), bad[:10])
# ON leaves and fold
Lv=leaves(); pins=json.load(open(os.path.join(HERE,'..','agentF_work','pins.json')))
on=[int(g) for g in pins if v[int(g)]%p==1]
print('ON leaf selectors:', sorted(on))
tot=None
for g in on: tot=add(tot,Lv[g])
print('fold =',((tot[0]-cs)%p, tot[1]))
print('L_24601 raw =',((Lv[24601][0]-cs)%p, Lv[24601][1]))
print('L_2081  raw =',((Lv[2081][0]-cs)%p, Lv[2081][1]))
print('x13682=',v[13682]%p,' x37892=',v[37892]%p)
print('x24468=',v[24468]%p,' x18956=',v[18956]%p)
