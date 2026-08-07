import sys, os, json, re, collections, pickle
F='/home/user/integer_solver/solve_lab/agentT_work/mirror/F'; sys.path.insert(0,F)
from fwd import Engine, NV
from parse import node_str
from circ2 import vars_of
E=Engine()
defrhs={c[1]:c[2] for c in (E.cls[a] for a in E.order)}
resby=collections.defaultdict(list)
for a in E.res:
    for u in vars_of(E.atoms[a]): resby[u].append(a)
uses=collections.defaultdict(list)
for w,r in defrhs.items():
    for u in vars_of(r): uses[u].append(w)
D=pickle.load(open('ortree2.pkl','rb'))
tree=D['tree']; selmap=D['selmap']

def isconst(v):
    r=defrhs.get(v)
    return r is not None and r[0]=='c'
def constval(v):
    r=defrhs.get(v)
    return r[1] if (r is not None and r[0]=='c') else None

leaves=[v for v in tree if tree[v] is None]
zero=[v for v in leaves if constval(v)==0]
free=[v for v in leaves if v not in defrhs]
other=[v for v in leaves if v not in zero and v not in free]
print('leaves %d : hardzero %d, free %d, other %d'%(len(leaves),len(zero),len(free),len(other)))
print('other sample',[(v,node_str(defrhs[v])[:60]) for v in other[:5]])

# selector -> gated products  sel*w
gated=collections.defaultdict(list)
for w,r in defrhs.items():
    if r[0]=='*' and r[1][0]=='v' and r[2][0]=='v':
        gated[r[1][1]].append((w,r[2][1])); gated[r[2][1]].append((w,r[1][1]))

stats=collections.Counter()
NODE={}
for n,ch in tree.items():
    if ch is None: continue
    sm=selmap[n]
    sa=sm['s_a']; sb=sm['s_b']; sab=sm['s_ab']
    # sab may contain the AND used inside the OR formula (used only in the OR def)
    sab_real=[s for s in sab if any(w!=n for w,_ in gated.get(s,[]))]
    ga=[g for s in sa for g in gated.get(s,[])]
    gb=[g for s in sb for g in gated.get(s,[])]
    gab=[g for s in sab_real for g in gated.get(s,[])]
    stats[(len(ga),len(gb),len(gab))]+=1
    NODE[n]=dict(a=sm['a'],b=sm['b'],sa=sa,sb=sb,sab=sab_real,ga=ga,gb=gb,gab=gab)
print('gated-term counts (a,b,ab):',stats.most_common(10))
pickle.dump(NODE,open('nodes.pkl','wb'))
