import os, sys, time
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentG_work')
import gsym
from gsym import *
src = sys.argv[1] if len(sys.argv)>1 else '/home/user/integer_solver/solve_lab/s10/AG_39013.json'
SYMS = [int(s) for s in sys.argv[2].split(',')] if len(sys.argv)>2 else [22162,30213]
v = L.load(src); ad.fwd(v, rounds=6)
av = L.all_atom_values(v); fail=L.failing_eqs(av)
print('base state %s score %d' % (src, L.NEQ-len(fail)))
val = gsym.build(v, SYMS, cap=int(sys.argv[3]) if len(sys.argv)>3 else 12)
# every check atom: symbolic polynomial
n=len(SYMS)
nz_sym=[]
for a in gsym.check_atoms():
    f = gsym.evalpoly_sym(a, val, n, 12)
    if isinstance(f,int):
        if f % P: nz_sym.append((a,'CONST_NONZERO',f))
    else:
        nz_sym.append((a,'SYM',f))
print('checks that are symbolic or nonzero-const:', len(nz_sym))
for a,kind,f in nz_sym[:60]:
    if kind=='SYM':
        print('  a%-6d deg=%d terms=%d  neq=%d' % (a,gsym.deg(f),len(f),len(L.atom2eq.get(a,{}))))
    else:
        print('  a%-6d CONST NONZERO neq=%d' % (a,len(L.atom2eq.get(a,{}))))
import pickle
pickle.dump({'syms':SYMS,'nz':nz_sym}, open('/home/user/integer_solver/solve_lab/agentG_work/sym_%s.pkl'%('_'.join(map(str,SYMS))[:40]),'wb'))
