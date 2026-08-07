import os, sys
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentG_work')
import gsym
from gsym import *
src = sys.argv[1] if len(sys.argv)>1 else '/home/user/integer_solver/solve_lab/s10/AG_39013.json'
SYMS = [int(s) for s in sys.argv[2].split(',')] if len(sys.argv)>2 else [22162,30213]
CAP = int(sys.argv[3]) if len(sys.argv)>3 else 12
v = L.load(src); ad.fwd(v, rounds=6)
print('score', L.NEQ-len(L.failing_eqs(L.all_atom_values(v))))
val = gsym.build(v, SYMS, cap=CAP)
n=len(SYMS)
def mstr(m):
    s=[]
    for i,e in enumerate(m):
        if e: s.append('s%d^%d'%(SYMS[i],e) if e>1 else 's%d'%SYMS[i])
    return '*'.join(s) if s else '1'
rows=[]
for a in gsym.check_atoms():
    f = gsym.evalpoly_sym(a, val, n, CAP)
    if isinstance(f,int):
        if f%P: print('a%-6d neq=%-3d CONST %d'%(a,len(L.atom2eq.get(a,{})),f%P))
        continue
    print('a%-6d neq=%-3d deg%d : %s'%(a,len(L.atom2eq.get(a,{})),gsym.deg(f),
        ' + '.join('%d*%s'%(c,mstr(m)) for m,c in sorted(f.items()))[:400]))
    rows.append((a,f))
