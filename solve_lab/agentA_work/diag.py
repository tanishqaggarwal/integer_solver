import sys, json, collections; sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentA_work')
from fractions import Fraction as F
import env, lib as L
from regsolve import analyse, qsolve
P=env.P
path=sys.argv[1]
v=L.load(path); av=L.all_atom_values(v); fe=L.failing_eqs(av)
nz=[a for a in range(L.NA) if av[a]]
print('%s score=%d failing=%d nz=%s'%(path.split('/')[-1],L.NEQ-len(fe),len(fe),nz))
A=sorted(set(a for e in fe for a in L.eq_atoms[e][2]))
K,R,rows,QUAD=analyse(v,A)
sol,free,incons,r,eqs=qsolve(rows,len(K))
bad=[(K[j],sol[j]) for j in range(len(K)) if sol[j] is not None and sol[j].denominator!=1]
print('atoms=%d knobs=%d eqs=%d rank=%d free=%d incons=%d'%(len(A),len(K),len(R),r,len(free),len(incons)))
for u,s in bad:
    print('  BAD knob x%d den=%s'%(u, 'p' if s.denominator==P else s.denominator))
    print('     atoms containing it:', [(a,'G%d'%L.atom_out[a][1] if a in L.atom_out else 'CHK', len(L.atom2eq[a]), 'NZ' if av[a] else '0') for a in L.var_atoms[u]])
    for a in L.var_atoms[u]:
        Pp=L.polys[a]
        if len(Pp)<=25:
            ts=[]
            for m,c in sorted(Pp.items()):
                cs='%+d'%c if abs(c)<10**12 else '%+d~%dd'%(c%1000,len(str(abs(c))))
                ts.append(cs+''.join('*x%d'%w for w in m))
            print('       a%d = %s'%(a,' '.join(ts)))
        else: print('       a%d [%d terms]'%(a,len(Pp)))
    print('     numerator mod p =', s.numerator % P, ' need p | numerator')
print('nonzero atoms and their required values:')
