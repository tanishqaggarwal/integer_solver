import sys, json, re, collections, itertools
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentC_work')
from close3 import *
from ort import leaves
P=2**256-2**32-977
W='/home/user/integer_solver/solve_lab/agentC_work/'
lp={int(k):[int(v[0]),int(v[1]),int(v[2]),int(v[3])] for k,v in json.load(open(W+'leafpts2.json')).items()}
side={}
for r,tag in [(8599,'s1'),(21839,'s1'),(25956,'s2'),(7304,'s2')]:
    for x in leaves(r):
        if x in lp: side[x]=tag
info={}
for b in sorted(lp):
    px,py,vx,vy=lp[b]
    pins={}
    for a in L.var_atoms[b]:
        if a in L.atom_out: continue
        m=re.match(r'^x_%d \* \(x_(\d+) - (-?\d+)\) - (?:(\d+) \* )?x_(\d+)$'%b, L.atom_src[a])
        if m:
            X=int(m.group(1))
            pins[X]=dict(atom=a,X=X,C=int(m.group(2)),mult=int(m.group(3) or 1),H=int(m.group(4)))
    if vx in pins and vy in pins:
        hx=outs.get(pins[vx]['H']); hy=outs.get(pins[vy]['H'])
        if hx is None or hy is None: continue
        info[b]=dict(side=side[b],px=px,py=py,
                     xpin=pins[vx],ypin=pins[vy],hx=hx,hy=hy,
                     Ex=set(L.atom2eq.get(hx,{})),Ey=set(L.atom2eq.get(hy,{})))
print('bits with full info:',len(info))
S1=[b for b in info if info[b]['side']=='s1']; S2=[b for b in info if info[b]['side']=='s2']
print('s1',len(S1),'s2',len(S2))
res=[]
for b,d in info.items():
    res.append((len(d['Ex']|d['Ey']),'A',b,b))
for u in S1:
    for w in S2:
        res.append((len(info[u]['Ex']|info[w]['Ey']),'Bxy',u,w))
        res.append((len(info[u]['Ey']|info[w]['Ex']),'Byx',u,w))
res.sort()
print('best 25 override plans (cost = |union of broken handle equations|):')
for r in res[:25]:
    print('   cost=%-3d %s  u=x_%-6d w=x_%-6d'%r)
print('single-handle eq counts: min Ex/Ey per side')
for tag,S in [('s1',S1),('s2',S2)]:
    a=sorted((len(info[b]['Ex']),b,'x') for b in S)[:5]
    c=sorted((len(info[b]['Ey']),b,'y') for b in S)[:5]
    print('  ',tag,'cheapest x-handles',a,' y-handles',c)
json.dump({str(b):{'side':d['side'],'hx':d['hx'],'hy':d['hy'],
                   'xpin':d['xpin'],'ypin':d['ypin'],'px':str(d['px']),'py':str(d['py'])} for b,d in info.items()},
          open(W+'bitinfo.json','w'))
