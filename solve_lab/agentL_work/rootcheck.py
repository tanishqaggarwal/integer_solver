import json, pickle, sys, collections, importlib.util
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentF_work')
from fwd import Engine,NV
E=Engine()
p=115792089237316195423570985008687907853269984665640564039457584007908834671663
M=pickle.load(open('full_model.pkl','rb')); OUT=M['OUT']; ROOT=M['ROOT']; NODE=M['NODE']
D=json.load(open('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json'))
vv=[0]*NV
for k,v in D.items():
    if k.startswith('x_'): vv[int(k[2:])]=int(v)
V=E.values(vv) if hasattr(E,'values') else None
def val(i):
    return vv[i]
# root mux input slot wires and selectors
d0,d1=OUT[ROOT]
print('ROOT x%d  va wires (%d,%d)  vb wires (%d,%d)  vab wires (%d,%d)'%(
    ROOT,d0['va'],d1['va'],d0['vb'],d1['vb'],d0['vab'],d1['vab']))
for nm,w in (('va0',d0['va']),('va1',d1['va']),('vb0',d0['vb']),('vb1',d1['vb']),
             ('vab0',d0['vab']),('vab1',d1['vab'])):
    print('   %-5s x%-6d = %s'%(nm,w,val(w)%p))
A=(val(d0['va'])%p,val(d1['va'])%p); B=(val(d0['vb'])%p,val(d1['vb'])%p)
print('A (root va input) == B (root vb input) ?', A==B)
print('  coord0 equal?',A[0]==B[0],'  coord1 equal?',A[1]==B[1])
sel=NODE[ROOT]['sab'][0]; sa=NODE[ROOT]['sa'][0]; sb=NODE[ROOT]['sb'][0]
print('root selectors: sel_a x%d=%s  sel_b x%d=%s  sel_ab x%d=%s'%(sa,val(sa),sb,val(sb),sel,val(sel)))
# what my model says the fold of {2081,24601} is
spec=importlib.util.spec_from_file_location('ss','/home/user/integer_solver/solve_lab/agentL_work/subsearch.py')
ss=importlib.util.module_from_spec(spec); spec.loader.exec_module(ss)
print('my fold2(2081,24601) =',ss.fold2(2081,24601))
print('TARGET               =',ss.TGT)
print('my model: leaf 24601 value at root frame =',ss.sw(ss.LEAF[24601],ss.cums[24601][-1]))
print('my model: leaf 2081  value at root frame =',ss.sw(ss.LEAF[2081],ss.cums[2081][-1]))
