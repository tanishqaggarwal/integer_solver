import sys; sys.path.insert(0,'.')
import env, lib as L
P=env.P
v=env.load_best(); av=L.all_atom_values(v)
def fmt(Pp,maxt=40):
    ts=[]
    for i,(m,c) in enumerate(sorted(Pp.items())):
        if i>=maxt: ts.append('...(%d more)'%(len(Pp)-maxt)); break
        cs = '%+d'%c if abs(c)<10**15 else '%+d[%dd]'%(c%1000, len(str(abs(c))))
        ts.append(cs + ''.join('*x%d'%u for u in m))
    return ' '.join(ts)
for a in [22231,37887,29090,39166,40066,40932,2202,16897,21113,38521,40005,40121,1465,8263,36088,1459,8261]:
    o = 'G%d'%L.atom_out[a][1] if a in L.atom_out else 'CHK'
    print('a%-6d [%s] val=%s eqs=%s'%(a,o,'NZ' if av[a] else '0', sorted(L.atom2eq[a])))
    print('     = %s'%fmt(L.polys[a]))
