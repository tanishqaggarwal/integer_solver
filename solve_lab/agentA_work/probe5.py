import sys, collections; sys.path.insert(0,'.')
import env, lib as L
P=env.P
v=env.load_best(); av=L.all_atom_values(v)
E=sorted(set(e for a in env.SEVEN for e in L.atom2eq[a])); Es=set(E)
A=sorted(set(a for e in E for a in L.eq_atoms[e][2]))
As=set(A)
V=sorted(set(u for a in A for u in L.avars[a]))
ext={u:[a for a in L.var_atoms[u] if a not in As] for u in V}
SAFE=[u for u in V if not ext[u]]
print('SAFE vars (0 external atoms):',SAFE)
def fmt(Pp,maxt=60):
    ts=[]
    for i,(m,c) in enumerate(sorted(Pp.items())):
        if i>=maxt: ts.append('...(%d more)'%(len(Pp)-maxt)); break
        cs='%+d'%c if abs(c)<10**12 else '%+d~%dd'%(c%1000,len(str(abs(c))))
        ts.append(cs+''.join('*x%d'%u for u in m))
    return ' '.join(ts)
for a in A:
    if len(L.polys[a])>25: 
        print('a%-6d [%d terms] vars=%s'%(a,len(L.polys[a]),sorted(L.avars[a])[:20]));continue
    print('a%-6d = %s'%(a,fmt(L.polys[a])))
