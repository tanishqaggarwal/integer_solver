"""Direct test of the cancellation hypothesis at the point of interest:
for each FAILING equation, list every atom in it with (#equations it occurs in,
private-handle granularity, current value).  An atom occurring in exactly ONE equation
and carrying a free private handle would make that equation trivially satisfiable."""
import sys, collections; sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentA_work')
import env, lib as L, ahandles as HH
P=env.P
v=L.load('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json')
av=L.all_atom_values(v); fe=L.failing_eqs(av)
solo,gran=HH.build(v)
one_eq=[a for a in range(L.NA) if len(L.atom2eq.get(a,{}))==1]
print('atoms occurring in exactly one equation: %d'%len(one_eq))
free1=[a for a in one_eq if gran.get(a)==1]
freep=[a for a in one_eq if gran.get(a)==P]
print('   of those, with a granularity-1 private handle: %d ; granularity-p: %d ; none: %d'%(
      len(free1),len(freep),len(one_eq)-len(free1)-len(freep)))
print()
for e in fe:
    m,sq,co=L.eq_atoms[e]
    print('eq%-7d mult=%-10d sq=%-5s atoms=%d'%(e,m,sq,len(co)))
    for a,c in sorted(co.items()):
        g=gran.get(a); gs='1' if g==1 else ('p' if g==P else ('-' if not g else str(g)))
        print('    a%-6d coeff%+4d  #eqs=%-3d handle_gran=%-3s value=%s'%(
            a,c,len(L.atom2eq[a]),gs,'NONZERO' if av[a] else '0'))
print()
# global: which equations COULD be fixed for free by a single-equation free atom?
fixable=set()
for a in free1:
    fixable |= set(L.atom2eq[a])
print('equations containing a single-equation granularity-1 atom: %d'%len(fixable))
print('   any of the 7 failing among them? %s'%sorted(set(fe)&fixable))
fixp=set()
for a in freep: fixp |= set(L.atom2eq[a])
print('equations containing a single-equation granularity-p atom: %d ; failing among them: %s'%(
      len(fixp),sorted(set(fe)&fixp)))
