import sys, collections; sys.path.insert(0,'.')
import env, lib as L
P=env.P
v=env.load_best(); av=L.all_atom_values(v)
E=sorted(set(e for a in env.SEVEN for e in L.atom2eq[a])); Es=set(E)
A=set()
for e in E: A|=set(L.eq_atoms[e][2])
# closure: equations touched by A
R=set()
for a in A: R|=set(L.atom2eq[a])
A2=set()
for e in R: A2|=set(L.eq_atoms[e][2])
print('level1 atoms %d -> eqs %d -> level2 atoms %d'%(len(A),len(R),len(A2)))
V=set()
for a in A: V|=L.avars[a]
print('variables in the 33 region atoms: %d'%len(V))
rows=[]
for u in sorted(V):
    outs=[a for a in L.var_atoms[u] if a not in A]
    rows.append((len(outs),u,outs))
rows.sort()
print('%-8s %-6s %-8s %s'%('var','#extA','value','outside atoms'))
for n,u,outs in rows:
    val=v[u]
    vs = str(val) if abs(val)<10**7 else ('p' if val==P else '%dd'%len(str(abs(val))))
    print('x%-7d %-6d %-8s %s'%(u,n,vs,outs[:12]))
