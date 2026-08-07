import sys, pickle, json
sys.path.insert(0,'/home/user/integer_solver/solve_lab')
sys.set_int_max_str_digits(2000000)
import checker
M=pickle.load(open('/home/user/integer_solver/solve_lab/agentT_work/mirror/L/full_model.pkl','rb'))
C=pickle.load(open('/home/user/integer_solver/solve_lab/agentT_work/mirror/L/calib2.pkl','rb'))
NODE=M['NODE']; OUT=M['OUT']; tree=M['tree']; live=set(M['live']); link=M['link']
sub=M['sub']; order=M['order']; PIN=M['PIN']; ROOT=M['ROOT']
ORIENT=C['ORIENT']; perm=C['perm']
p=115792089237316195423570985008687907853269984665640564039457584007908834671663
v0=checker.load_assignment('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json')
XY=pickle.load(open('w_xy.pkl','rb'))
print('PIN[24601]',PIN[24601][0], [x%10**15 for x in PIN[24601][1]])
print('deliv x_33462 mod p =', v0[33462]%p%10**15, ' x_22152 mod p =',v0[22152]%p%10**15)
print('X const(72) mod p =',XY['X'][72][3]%p%10**15,'  Y const(72) mod p=',XY['Y'][72][3]%p%10**15)
# how many nodes have nonzero vab in the deliverable?
nz=[]
for n in order:
    vals=[v0[d['vab']] for d in OUT[n]]
    if any(vals): nz.append(n)
print('nodes with nonzero vab in deliverable: %d'%len(nz))
print('ROOT',ROOT,'in nz?',ROOT in nz)
# count nonzero va/vb
cnt=0
for n in order:
    for d in OUT[n]:
        for s in ('va','vb'):
            if v0[d[s]]: cnt+=1
print('nonzero va/vb wires:',cnt,'of',383*4)
# depth of each node
depth={}
def dep(n):
    if n in depth: return depth[n]
    if tree.get(n) is None: depth[n]=0; return 0
    a,b=tree[n]; d=1+max(dep(a),dep(b)); depth[n]=d; return d
for n in order: dep(n)
print('ROOT depth',depth[ROOT])
import collections
print('live-slot nodes (both subtrees contain live leaves):',sum(1 for n in order if [x for x in sub[tree[n][0]] if x in live] and [x for x in sub[tree[n][1]] if x in live]))
pickle.dump({'depth':depth},open('x_depth.pkl','wb'))
