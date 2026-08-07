import sys; sys.path.insert(0,'.')
import env, lib as L
P=env.P
v=env.load_best(); av=L.all_atom_values(v)
SEVEN=env.SEVEN
E=sorted(set(e for a in SEVEN for e in L.atom2eq[a]))
Es=set(E)
VARS=[642,7068,2099,28730,9413,1329,29854,9118,7075,31864,10903,8731,17325,
      17499,22665,28961,28599]
print('values:')
for u in VARS:
    val=v[u]
    print('  x%-6d = %s%s'%(u, str(val)[:40]+('...' if len(str(val))>40 else ''),
          '   (== p)' if val==P else ('  (=%d)'%val if abs(val)<10**6 else '')))
print()
print('atom membership (atoms containing each var), with status:')
for u in VARS:
    ats=L.var_atoms[u]
    out=[]
    for a in ats:
        st='NZ' if av[a] else 'z'
        g='G%d'%L.atom_out[a][1] if a in L.atom_out else 'CHK'
        neq=len(L.atom2eq.get(a,{}))
        outside=len(set(L.atom2eq.get(a,{}))-Es)
        out.append('%d[%s,%s,eq%d,out%d]'%(a,st,g,neq,outside))
    print('  x%-6d (%d atoms): %s'%(u,len(ats),' '.join(out)))
