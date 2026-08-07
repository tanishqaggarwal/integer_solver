import sys, os, collections
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath('.')), 'solve_lab','s9','eff'))
sys.path.insert(0,'/home/user/integer_solver/solve_lab/s9/eff')
import lib as L
print('NA(atoms)=',L.NA,'NEQ=',L.NEQ,'NVARS=',L.NVARS)
print('gate atoms (have output):',len(L.atom_out))
chk=[a for a in range(L.NA) if a not in L.atom_out]
print('check atoms:',len(chk))
# free inputs: vars that are not any gate's output
outs=set(t for (c,t) in L.atom_out.values())
print('gate-output vars:',len(outs))
allv=set()
for s in L.avars: allv|=s
print('vars used:',len(allv))
free=allv-outs
print('free inputs:',len(free))
# atoms appearing in equations
inaeq=set(L.atom2eq)
print('atoms in equations:',len(inaeq))
print('check atoms in equations:',len(set(chk)&inaeq))
print('gate atoms in equations:',len(set(L.atom_out)&inaeq))
# degree distribution of atoms
dd=collections.Counter()
for Pp in L.polys:
    dd[max(len(m) for m in Pp)]+=1
print('atom degrees:',dict(dd))
# equation shapes
sq=sum(1 for (m,s,c) in L.eq_atoms if s)
print('square eqs:',sq,'linear eqs:',L.NEQ-sq)
sz=collections.Counter(len(c) for (m,s,c) in L.eq_atoms)
print('eq atom-count dist:',sorted(sz.items())[:15])
