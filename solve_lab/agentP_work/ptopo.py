#!/usr/bin/env python3
"""Agent P: definitive topology of the 382-stage network."""
import pickle,sys,json
from collections import Counter,defaultdict,deque
sys.set_int_max_str_digits(10**7)
W='/home/user/integer_solver/solve_lab/agentP_work/'
P=115792089237316195423570985008687907853269984665640564039457584007908834671663
D=pickle.load(open(W+'model4.pkl','rb')); AP=D['AP']
S=pickle.load(open(W+'slp.pkl','rb')); topo=S['topo']; outof=S['outof']
G=pickle.load(open(W+'graph.pkl','rb')); stages=G['stages']; outmap=G['outmap']
Wd=pickle.load(open(W+'wire3.pkl','rb')); par=Wd['par']; leafcoord=Wd['leafcoord']; Z=set(Wd['Z'])
R=lambda x:par[x]
# zero-pinned classes
zero=set()
for ap in AP:
    if len(ap)==1:
        m=next(iter(ap))
        if len(m)==1: zero.add(R(m[0]))
print("zero-pinned classes:",len(zero))

src=[]   # per stage: [(kindX,..),(kindY,..)]
for s in stages:
    row=[]
    for pr in (s['X'],s['Y']):
        a,b=R(pr[0]),R(pr[1])
        ka,kb=outmap.get(a),outmap.get(b)
        if ka and kb and ka[0]==kb[0]: row.append(('S',ka[0]))
        elif a in leafcoord and b in leafcoord and leafcoord[a][0]==leafcoord[b][0]: row.append(('L',leafcoord[a][0]))
        elif a in zero and b in zero: row.append(('0',None))
        else: row.append(('?',pr))
    src.append(row)
print("kinds:",Counter(k for r in src for k,_ in r))

consumed=Counter()
for r in src:
    for k,v in r:
        if k=='S': consumed[v]+=1
print("stage out-degree hist:",sorted(Counter(consumed[j] for j in range(len(stages))).items()))
roots=[j for j in range(len(stages)) if consumed[j]==0]
print("stages whose output feeds no stage:",len(roots),roots[:40])

# depth from leaves
depth=[None]*len(stages)
def dep(j,seen=None):
    if depth[j] is not None: return depth[j]
    d=0
    for k,v in src[j]:
        if k=='S': d=max(d,dep(v)+1)
        elif k=='L': d=max(d,1)
    depth[j]=d; return d
sys.setrecursionlimit(10000)
for j in range(len(stages)): dep(j)
print("depth hist:",sorted(Counter(depth).items()))

# leaf support per stage
supp=[None]*len(stages)
def sup(j):
    if supp[j] is not None: return supp[j]
    s=set()
    for k,v in src[j]:
        if k=='S': s|=sup(v)
        elif k=='L': s.add(v)
    supp[j]=s; return s
for j in range(len(stages)): sup(j)
print("leaf-support size hist:",sorted(Counter(len(supp[j]) for j in range(len(stages))).items())[:40])
print("max support:",max(len(x) for x in supp))
# which stage covers all 256?
full=[j for j in range(len(stages)) if len(supp[j])==256]
print("stages with full 256 support:",full)
pickle.dump({'src':src,'depth':depth,'supp':[sorted(x) for x in supp],'roots':roots},open(W+'topo.pkl','wb'))
