import os, sys, pickle, collections
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentG_work')
import gsym2 as G
from gsym2 import L, ad, P
D=pickle.load(open('/home/user/integer_solver/solve_lab/agentG_work/coset_model.pkl','rb'))
Lin=pickle.load(open('/home/user/integer_solver/solve_lab/agentG_work/coset_lin.pkl','rb'))
NB=D['NB']; pt=D['pt']; x0=Lin['x0']; col=Lin['col']; n=len(NB)
supp=[c for c in range(n) if (pt[c]-x0[c])%P]
print('deliverable departure support (%d unknowns):'%len(supp))
for c in supp:
    print('   x%-6d  occurs in %d linear equations %s'%(NB[c],len(col.get(c,[])),col.get(c,[])[:14]))
occ=[c for c in supp if col.get(c)]
print('\nof those, %d occur in linear equations at all; %d are free of them'%(len(occ),len(supp)-len(occ)))
U=set()
for c in occ: U|=set(col[c])
print('union of their linear footprints: %d equations'%len(U))
print(sorted(U))
