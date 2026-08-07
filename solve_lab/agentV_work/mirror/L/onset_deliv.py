"""Settle the deliverable's ON-set directly from its own assignment file."""
import json, pickle, sys, collections
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentV_work/mirror/F')
from fwd import Engine, NV
E=Engine()
p=115792089237316195423570985008687907853269984665640564039457584007908834671663
M=pickle.load(open('full_model.pkl','rb'))
live=M['live']; dead=M['dead']; tree=M['tree']; NODE=M['NODE']; ROOT=M['ROOT']; sub=M['sub']
D=json.load(open('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json'))
def get(i):
    v=D.get('x_%d'%i,0)
    return int(v) if not isinstance(v,int) else v
print('deliverable json: %d entries'%len(D))
# TEST 1: direct read of the 256 free leaf variables
on=[L for L in live if get(L)!=0]
print('TEST 1  free leaf vars set nonzero in the file: %s'%sorted(on))
for L in sorted(on): print('        x%-6d = %r'%(L,get(L)))
print('        (all 256 leaves are FREE vars, so the file value IS the bit; unset => 0)')
# TEST 2: propagate and read every OR node's live flag
vv=[0]*NV
for k,v in D.items():
    if k.startswith('x_'): vv[int(k[2:])]=int(v)
r=E.run(vv)
defv={c[1]:c[2] for c in (E.cls[a] for a in E.order)}
full=E.eval_all(vv) if hasattr(E,'eval_all') else None
print('TEST 2  is 2081 one of my 256 live leaves? %s   one of the 128 dead leaves? %s'%(2081 in set(live),2081 in set(dead)))
print('        is 2081 a free variable at all? %s'%(2081 not in defv))
print('        deliverable value of x2081 = %r'%get(2081))
print('        deliverable value of x24601 = %r'%get(24601))
# where does 2081 sit
if 2081 in set(live):
    n=[k for k in NODE if NODE[k]['a']==2081 or NODE[k]['b']==2081]
    print('        2081 is a leaf under node(s)',n)
else:
    ln={}
    for k in NODE:
        for s,ch in (('va',NODE[k]['a']),('vb',NODE[k]['b'])): ln.setdefault(ch,[]).append((k,s))
    print('        2081 as an OR-tree child:',ln.get(2081))
    print('        2081 in defrhs? ',2081 in defv, defv.get(2081))
