#!/usr/bin/env python3
"""AUDIT T25 -- agent O's eq8680 Lemma: 'eq8680 = T^2 with T a linear form in 20 atoms, so a
square has a single zero locus and every satisfying assignment has T = 0.'
Four results hang on it.  Check the factorisation directly against the raw file."""
import re,sys,random,os
LINE=8681
lhs=open('/home/user/integer_solver/EQUATIONS.txt').read().split('\n')[LINE-1].rsplit('=',1)[0].strip()
print('raw LHS length: %d chars'%len(lhs))
print('occurrences of x_4432 : %d'%len(re.findall(r'x_4432\b',lhs)))
print('occurrences of x_28730: %d'%len(re.findall(r'x_28730\b',lhs)))
print('occurrences of x_19964: %d'%len(re.findall(r'x_19964\b',lhs)))
def strip_outer(s):
    while s.startswith('(') and s.endswith(')'):
        d=0; ok=True
        for k,ch in enumerate(s):
            if ch=='(': d+=1
            elif ch==')':
                d-=1
                if d==0 and k!=len(s)-1: ok=False; break
        if not ok: break
        s=s[1:-1]
    return s
def split_top(s,op='*'):
    s=strip_outer(s); parts=[]; d=0; cur=''
    i=0
    while i<len(s):
        ch=s[i]
        if ch=='(': d+=1
        elif ch==')': d-=1
        if d==0 and ch==op:
            parts.append(cur); cur=''; i+=1; continue
        cur+=ch; i+=1
    parts.append(cur)
    return [p for p in parts]
facs=split_top(lhs,'*')
print('\ntop-level * factors: %d'%len(facs))
for i,f in enumerate(facs):
    print('   factor %d: len %-6d  starts %s'%(i,len(f),f[:46].replace(' ','')))
norm=[strip_outer(f).replace(' ','') for f in facs]
print('\nall factors textually identical? %s'%(len(set(norm))==1))
if len(set(norm))>1:
    from collections import Counter
    c=Counter(norm)
    print('   distinct factor texts: %d ; multiplicities %s'%(len(c),sorted(c.values(),reverse=True)))
# numeric identity test
NV=38748
VAR=re.compile(r'x_(\d+)')
code_lhs=compile(VAR.sub(r'v[\1]',lhs),'<lhs>','eval')
Ttxt=strip_outer(facs[0])
code_T=compile(VAR.sub(r'v[\1]',Ttxt),'<T>','eval')
rnd=random.Random(5)
print('\nnumeric test at random points: is LHS == T^k ?')
for trial in range(4):
    v=[0]*NV
    for u in {int(m) for m in VAR.findall(lhs)}: v[u]=rnd.randint(-40,40)
    ns={'v':v,'__builtins__':{}}
    L=eval(code_lhs,ns); Tv=eval(code_T,ns)
    ks=[k for k in range(1,7) if Tv**k==L]
    print('   trial %d: T=%-14s LHS matches T^k for k in %s'%(trial,str(Tv)[:14],ks))
# is T linear?  test additivity in each variable
vars_=sorted({int(m) for m in VAR.findall(Ttxt)})
print('\nT: %d distinct variables'%len(vars_))
v0=[0]*NV
for u in vars_: v0[u]=rnd.randint(-30,30)
base=eval(code_T,{'v':v0,'__builtins__':{}})
nonlin=[]
for u in vars_:
    ys=[]
    for t in (0,1,2):
        vv=list(v0); vv[u]=v0[u]+t
        ys.append(eval(code_T,{'v':vv,'__builtins__':{}}))
    if ys[2]-ys[1]!=ys[1]-ys[0]: nonlin.append(u)
print('   variables in which T is NOT affine: %d %s'%(len(nonlin),nonlin[:8]))
print('   dT/dx_4432  = %s'%(eval(code_T,{'v':[x+ (1 if i==4432 else 0) for i,x in enumerate(v0)],'__builtins__':{}})-base))
print('   dT/dx_28730 = %s'%(eval(code_T,{'v':[x+ (1 if i==28730 else 0) for i,x in enumerate(v0)],'__builtins__':{}})-base))
print('   dT/dx_19964 = %s'%(eval(code_T,{'v':[x+ (1 if i==19964 else 0) for i,x in enumerate(v0)],'__builtins__':{}})-base))
