import os, sys, json, pickle
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentG_work')
import gsym
from gsym import *
d=pickle.load(open('sys112.pkl','rb')); SYMS=d['syms']; rows=d['rows']
D=dict(rows); n=len(SYMS); ix={u:i for i,u in enumerate(SYMS)}
v=L.load('/home/user/integer_solver/solve_lab/s10/AG_39013.json'); ad.fwd(v,rounds=6)
# extract A' and B from a19297 (8646263 A + 1073965 B) and a19299 (10159099 A + 6926539 B)
f1,f2=D[19297],D[19299]
# build A and B polynomials from the pencil
c11,c12,c21,c22=8646263,1073965,10159099,6926539
det=(c11*c22-c12*c21)%P; inv=pow(det,-1,P)
def lc(f,g,a,b):
    out=0
    for m,c in f.items(): out=gsym.padd(out,{m:c*a%P},n)
    for m,c in g.items(): out=gsym.padd(out,{m:c*b%P},n)
    return {m:c%P for m,c in out.items() if c%P}
A=lc(f1,f2,c22*inv%P,(-c12)*inv%P)
B=lc(f1,f2,(-c21)*inv%P,c11*inv%P)
def mstr(m): return '*'.join(('x%d'%SYMS[i]+('^%d'%e if e>1 else '')) for i,e in enumerate(m) if e) or '1'
print('A (%d terms):'%len(A))
for m,c in sorted(A.items(),key=lambda kv:-sum(kv[0])): print('   %-32s %d'%(mstr(m),c))
print('B (%d terms):'%len(B))
for m,c in sorted(B.items(),key=lambda kv:-sum(kv[0])): print('   %-32s %d'%(mstr(m),c))
# pinned coordinate values from the linear checks
C={}
lin={a:f for a,f in rows if gsym.deg(f)==1}
def pinval(a):
    f=lin[a]; const=f.get((0,)*n,0); terms={m:c for m,c in f.items() if sum(m)}
    return const,terms
for a in [1618,688,3576,3578,31670,31672]:
    print('pin a%d ->'%a, pinval(a)[0], [(SYMS[[i for i,e in enumerate(m) if e][0]],c) for m,c in pinval(a)[1].items()])
