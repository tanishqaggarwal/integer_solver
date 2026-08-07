import sys, collections; sys.path.insert(0,'.')
import env, lib as L, ahandles as HH
P=env.P
v=env.load_best(); av=L.all_atom_values(v)
solo,gran=HH.build(v)
S=set(env.SEVEN)
E=sorted(set(e for a in env.SEVEN for e in L.atom2eq[a])); Es=set(E)
A=set()
for e in E: A|=set(L.eq_atoms[e][2])
print('atoms in region E: %d'%len(A))
print('%-8s %-6s %-5s %-6s %-6s %s'%('atom','role','#eq','outE','gran','out-equations'))
for a in sorted(A):
    out=sorted(set(L.atom2eq[a])-Es)
    role='S' if a in S else ('G%d'%L.atom_out[a][1] if a in L.atom_out else 'CHK')
    g=gran.get(a); gs='p' if g==P else (str(g) if g else '-')
    print('%-8d %-6s %-5d %-6d %-6s %s'%(a,role,len(L.atom2eq[a]),len(out),gs,out[:10]))
print()
m,sq,co=L.eq_atoms[8680]
print('eq8680 mult=%d sq=%s atoms:'%(m,sq))
for a,c in sorted(co.items()):
    g=gran.get(a); gs='p' if g==P else (str(g) if g else '-')
    role='G%d'%L.atom_out[a][1] if a in L.atom_out else 'CHK'
    print('   a%-6d %+d  role=%s neq=%d gran=%s val=%s'%(a,c,role,len(L.atom2eq[a]),gs,'NZ' if av[a] else '0'))
