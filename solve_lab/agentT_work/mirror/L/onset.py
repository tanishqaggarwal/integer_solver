import sys, os, json, re, collections, pickle
F='/home/user/integer_solver/solve_lab/agentT_work/mirror/F'; sys.path.insert(0,F)
from fwd import Engine, NV
from parse import node_str
from circ2 import vars_of
E=Engine()
p=115792089237316195423570985008687907853269984665640564039457584007908834671663
defrhs={c[1]:c[2] for c in (E.cls[a] for a in E.order)}
D=pickle.load(open('ortree2.pkl','rb')); tree=D['tree']
NODE=pickle.load(open('nodes.pkl','rb')); OUT=pickle.load(open('outwires.pkl','rb'))
leaves=[v for v in tree if tree[v] is None]
freeleaf=[v for v in leaves if v not in defrhs]
A=json.load(open('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json'))
v=[0]*NV
for k,x in A.items(): v[int(k[2:])]=int(x)
E.run(v)  # propagate defined vars
on=[l for l in freeleaf if v[l]!=0]
print('free leaves',len(freeleaf),' ON in deliverable:',len(on), on)
print('values of ON leaves',[v[l] for l in on])
# which nodes are live
def live(n):
    if tree[n] is None: return v[n]%p!=0
    return live(tree[n][0]) or live(tree[n][1])
print('nodes live count', sum(1 for n in NODE if live(n)))
pickle.dump({'freeleaf':freeleaf,'on':on},open('leafinfo.pkl','wb'))
