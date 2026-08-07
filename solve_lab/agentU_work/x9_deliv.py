import sys, pickle
sys.path.insert(0,'/home/user/integer_solver/solve_lab')
sys.set_int_max_str_digits(2000000)
import checker
M=pickle.load(open('/home/user/integer_solver/solve_lab/agentT_work/mirror/L/full_model.pkl','rb'))
C=pickle.load(open('/home/user/integer_solver/solve_lab/agentT_work/mirror/L/calib2.pkl','rb'))
NODE=M['NODE']; OUT=M['OUT']; tree=M['tree']; live=set(M['live']); link=M['link']
sub=M['sub']; order=M['order']; PIN=M['PIN']; ROOT=M['ROOT']; leafnode=M['leafnode']
v0=checker.load_assignment('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json')
L=pickle.load(open('v_leaves.pkl','rb')); sel2exp=L['sel2exp']
D=pickle.load(open('x_depth.pkl','rb'))['depth']
# path from leaf-node to root
par={}
for n in tree:
    if tree[n]: 
        for ch in tree[n]: par[ch]=n
def path(x):
    q=[]; 
    while x in par: q.append(par[x]); x=par[x]
    return q
print('== nodes with nonzero vab ==')
for n in tree:
    if tree[n] is None: continue
    vals=[v0[d['vab']] for d in OUT[n]]
    if any(vals):
        a,b=tree[n]
        print(' node %d depth=%d |sub|=%d  vab bits=%s'%(n,D[n],len(sub[n]),[x.bit_length() for x in vals]))
        print('    va=%s vb=%s'%([v0[d['va']].bit_length() for d in OUT[n]],[v0[d['vb']].bit_length() for d in OUT[n]]))
        print('    vab vals:', [str(x)[:50] for x in vals])
print()
print('== leaf nodes of the two ON selectors ==')
for s in (2081,24601):
    ws=PIN[s][0]
    ln=[k for k,vv in leafnode.items() if k in ws]
    print(' sel %d exp %d wires %s'%(s,sel2exp[s],ws))
    for w in ws:
        if w in leafnode: print('    wire %d -> leafnode %s ; path %s'%(w,leafnode[w],path(leafnode[w][0])[:12]))
print()
print('== ROOT out wires and values ==')
for i,d in enumerate(OUT[ROOT]):
    print('  i=%d va=%d(%d bits) vb=%d(%d bits) vab=%d(%d bits) out=%d(%d bits)'%(
      i,d['va'],v0[d['va']].bit_length(),d['vb'],v0[d['vb']].bit_length(),
      d['vab'],v0[d['vab']].bit_length(),d['out'],v0[d['out']].bit_length()))
