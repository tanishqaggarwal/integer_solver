import pickle, sys, collections
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentV_work/mirror/F')
from fwd import Engine,NV
E=Engine()
M=pickle.load(open('full_model.pkl','rb'))
eqs_of=collections.defaultdict(set)
for e,row in enumerate(E.eqres):
    for k,j in row: eqs_of[E.res[j]].add(e)
cnt=collections.Counter({a:len(s) for a,s in eqs_of.items()})
print('atom eq-count hist',sorted(collections.Counter(cnt.values()).items())[:12])
link=M['link']; OUT=M['OUT']; NODE=M['NODE']; tree=M['tree']
# atom text for a slot wire w
atom_of={}
for a in E.res:
    pass
import re
lr=[re.compile(r'^\(\(x(\d+)-x(\d+)\)[-+]'), re.compile(r'^\(\((\d+)\*\(x(\d+)-x(\d+)\)\)-')]
for a in E.res:
    m=lr[0].match(a)
    if m: u,z=int(m.group(1)),int(m.group(2))
    else:
        m=lr[1].match(a)
        if not m: continue
        u,z=int(m.group(2)),int(m.group(3))
    for w in (u,z):
        if w in link: atom_of[w]=a
rows=[]
for n in NODE:
    for side,ch in (('va',NODE[n]['a']),('vb',NODE[n]['b'])):
        if tree[ch] is None: continue
        ws=[d[side] for d in OUT[n]]
        if not all(w in atom_of for w in ws): continue
        u=set()
        for w in ws: u|=eqs_of[atom_of[w]]
        rows.append((len(u),n,side,ws))
rows.sort()
print('slot-link pairs priced:',len(rows))
print('cheapest 15:')
for r in rows[:15]: print('   cost %-3d node x%-6d %s wires %s'%r)
print('cost histogram',sorted(collections.Counter(r[0] for r in rows).items())[:12])
# also the top/root pair used by mkassign
print('top pair (24468,18956) cost', len(eqs_of[atom_of[24468]]|eqs_of[atom_of[18956]]) if 24468 in atom_of else 'n/a')
pickle.dump({'rows':rows,'atom_of':atom_of,'eqs_of':{a:sorted(s) for a,s in eqs_of.items()}},open('price.pkl','wb'))
