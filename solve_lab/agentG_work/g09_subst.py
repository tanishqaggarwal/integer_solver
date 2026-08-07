import os, sys, json, pickle
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentG_work')
import gsym
from gsym import *
d=pickle.load(open('linpar.pkl','rb'))
M,piv,free,SYMS,non=d['M'],d['piv'],d['free'],d['SYMS'],d['nonlin']
n=len(SYMS); k=len(free)
# substitution: each original var index c -> poly in the k free params
sub=[None]*n
for j,c in enumerate(free):
    e=[0]*k; e[j]=1
    sub[c]={tuple(e):1}
for r,c in enumerate(piv):
    f={}
    const=M[r][n]%P
    if const: f[(0,)*k]=const
    for j,c2 in enumerate(free):
        co=(-M[r][c2])%P
        if co:
            e=[0]*k; e[j]=1
            f[tuple(e)]=co
    sub[c]=f if f else 0
def subpoly(f):
    out=0
    for m,c in f.items():
        t=c%P
        for i,e in enumerate(m):
            for _ in range(e):
                t=gsym.pmul(t,sub[i],k,None)
                if t==0: break
            if t==0: break
        if t!=0: out=gsym.padd(out,t,k)
    return out
res=[]
for a,f in non:
    g=subpoly(f)
    res.append((a,g))
    print('a%-6d -> deg %s terms %s' % (a, gsym.deg(g), gsym.nterms(g)))
pickle.dump({'free':[SYMS[c] for c in free],'polys':res,'sub':sub,'SYMS':SYMS,'piv':piv,'M':M,'freeidx':free}, open('redsys.pkl','wb'))
# print them
names=['t%d'%SYMS[c] for c in free]
def mstr(m):
    s=[]
    for i,e in enumerate(m):
        if e: s.append(names[i]+('^%d'%e if e>1 else ''))
    return '*'.join(s) if s else '1'
for a,g in res:
    if isinstance(g,int):
        print('a%d = CONST %d'%(a,g)); continue
    print('\na%d:'%a)
    for m,c in sorted(g.items(), key=lambda kv:(-sum(kv[0]),kv[0])):
        print('   %-30s %d'%(mstr(m),c))
