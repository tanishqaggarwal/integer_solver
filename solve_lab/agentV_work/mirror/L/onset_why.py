import json, pickle, sys, collections
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentT_work/mirror/F')
from fwd import Engine,NV
E=Engine()
M=pickle.load(open('full_model.pkl','rb')); NODE=M['NODE']; ROOT=M['ROOT']; sub=M['sub']; tree=M['tree']
live=set(M['live'])
try:
    O2=pickle.load(open('ortree2.pkl','rb'))
    pl=set(O2['live']) if isinstance(O2,dict) and 'live' in O2 else None
    print('stale partial model live-leaf count:',len(pl) if pl else '?')
    print('  2081 in stale partial live set?',2081 in pl if pl else '?')
    print('  24601 in stale partial live set?',24601 in pl if pl else '?')
except Exception as e: print('stale model:',e)
parent={}
for n in NODE:
    for s,ch in (('va',NODE[n]['a']),('vb',NODE[n]['b'])): parent[ch]=n
def anc(x):
    r=[]; 
    while x!=ROOT: x=parent[x]; r.append(x)
    return r
a1=anc(2081); a2=anc(24601)
s1=set(a1)
L=[n for n in a2 if n in s1][0]
print('LCA(2081,24601) = x%d ; is it the ROOT? %s'%(L,L==ROOT))
print('depth of 2081 path %d, of 24601 path %d'%(len(a1),len(a2)))
ra,rb=tree[ROOT]
print('root children x%d x%d ; 2081 under %s ; 24601 under %s'%(ra,rb,
   'a' if 2081 in sub[ra] else 'b', 'a' if 24601 in sub[ra] else 'b'))
