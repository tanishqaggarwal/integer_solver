import os, sys, json, pickle
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentG_work')
import gsym
from gsym import *
d=pickle.load(open('sys112.pkl','rb')); SYMS=d['syms']; rows=d['rows']
def mstr(m):
    s=[]
    for i,e in enumerate(m):
        if e: s.append('x%d'%SYMS[i]+('^%d'%e if e>1 else ''))
    return '*'.join(s) if s else '1'
def sc(c):
    c%=P
    return str(c) if c < 10**12 else ('(-%d)'%(P-c) if P-c < 10**12 else 'C[%s..]'%str(c)[:8])
for a,f in rows:
    if gsym.deg(f)>1:
        print('\n=== a%d (deg %d, %d eqs)'%(a,gsym.deg(f),len(L.atom2eq.get(a,{}))))
        for m,c in sorted(f.items(), key=lambda kv:(-sum(kv[0]),kv[0])):
            print('   %-40s %s'%(mstr(m),sc(c)))
print('\n=== LINEAR checks (%d) ==='%sum(1 for a,f in rows if gsym.deg(f)==1))
for a,f in rows:
    if gsym.deg(f)==1:
        terms=' '.join('%s*%s'%(sc(c),mstr(m)) for m,c in sorted(f.items(),key=lambda kv:-sum(kv[0])))
        print('a%-6d(%2d eq): %s'%(a,len(L.atom2eq.get(a,{})),terms[:170]))
