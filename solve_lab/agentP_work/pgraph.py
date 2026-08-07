#!/usr/bin/env python3
"""Agent P: definitive stage graph. stage j = law at q_j, output = mux at q_{j-1}+28.."""
import pickle,sys,json
from collections import Counter,defaultdict,deque
sys.set_int_max_str_digits(10**7)
W='/home/user/integer_solver/solve_lab/agentP_work/'
D=pickle.load(open(W+'model4.pkl','rb')); AP=D['AP']
S=pickle.load(open(W+'slp.pkl','rb')); topo=S['topo']; outof=S['outof']
B=pickle.load(open(W+'blocks.pkl','rb'))
Wd=pickle.load(open(W+'wire3.pkl','rb')); par=Wd['par']; leafcoord=Wd['leafcoord']
qpos=pickle.load(open(W+'qpos.pkl','rb'))['qpos']
R=lambda x: par[x]

def prod(ap):
    o=None;k=None;uv=None
    for m,c in ap.items():
        if len(m)==1 and o is None and abs(c)==1: o=(m[0],c)
        elif len(m)==2: uv=m; k=c
        else: return None
    if o is None or uv is None: return None
    return (o[0], k, uv[0], uv[1])

# stage j: law block index j (0..381 matched), output mux at q_{j}-43+28
stages=[]
for j,b in enumerate(B):
    mq=b['q']-43+28
    try:
        o1=outof[topo[mq+4]]; o2=outof[topo[mq+9]]
        s0=[prod(AP[topo[mq+t]]) for t in range(3)]
        s1=[prod(AP[topo[mq+5+t]]) for t in range(3)]
        ok = ({R(s0[0][2]),R(s0[0][3])} >= {R(b['i2'])} and
              {R(s0[1][2]),R(s0[1][3])} >= {R(b['i1'])} and
              {R(s0[2][2]),R(s0[2][3])} >= {R(b['i5'])} and
              {R(s1[0][2]),R(s1[0][3])} >= {R(b['i3'])} and
              {R(s1[1][2]),R(s1[1][3])} >= {R(b['i4'])} and
              {R(s1[2][2]),R(s1[2][3])} >= {R(b['i6'])})
    except Exception as e:
        ok=False; o1=o2=None
    stages.append(dict(j=j,q=b['q'],X=(b['i2'],b['i3']),Y=(b['i1'],b['i4']),Zlaw=(b['i5'],b['i6']),
                       out=(o1,o2),ok=ok,live=b['live']))
print("stages with a matching mux:",sum(1 for s in stages if s['ok']),"of",len(stages))

outmap={}
for s in stages:
    if s['out'][0] is not None:
        outmap[R(s['out'][0])]=(s['j'],0); outmap[R(s['out'][1])]=(s['j'],1)
print("stage output classes:",len(outmap))
kinds=Counter(); parents=defaultdict(list)
for s in stages:
    for nm,pr in (('X',s['X']),('Y',s['Y'])):
        a,b2=R(pr[0]),R(pr[1])
        ka = outmap.get(a); kb=outmap.get(b2)
        if ka is not None and kb is not None and ka[0]==kb[0]:
            kinds['stage']+=1; parents[s['j']].append(('S',ka[0]))
        elif a in leafcoord and b2 in leafcoord and leafcoord[a][0]==leafcoord[b2][0]:
            kinds['leaf']+=1; parents[s['j']].append(('L',leafcoord[a][0]))
        else:
            kinds['unknown']+=1; parents[s['j']].append(('?',pr))
print("stage input pair kinds:",kinds)
pickle.dump({'stages':stages,'parents':dict(parents),'outmap':outmap},open(W+'graph.pkl','wb'))
