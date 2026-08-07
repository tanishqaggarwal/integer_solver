#!/usr/bin/env python3
"""Spine-aware atom decomposition: S = A0 + c1*A1 + c2*A2 + ..."""
import re, sys, pickle, time, os, collections
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE)
from parse import parse_line, const_val, node_str
from parse2 import core_of, factors
EQ = os.path.join(HERE, '..', '..', 'EQUATIONS.txt')

def split_coef(nd):
    fs=[]; factors(nd,fs); k=1; nonc=[]
    for f in fs:
        cv=const_val(f)
        if cv is None: nonc.append(f)
        else: k*=cv
    if not nonc: return k, ('c',1)
    if len(nonc)==1: return k, nonc[0]
    n2=nonc[0]
    for f in nonc[1:]: n2=('*',n2,f)
    return k,n2

def spine(S):
    terms=[]; n=S
    while n[0] in ('+','-'):
        sg = 1 if n[0]=='+' else -1
        k,a = split_coef(n[2])
        terms.append((sg*k,a))
        n=n[1]
    k,a=split_coef(n)
    terms.append((k,a))
    terms.reverse()
    return terms

if __name__=='__main__':
    t0=time.time()
    lines=open(EQ).read().splitlines()
    eqs=[]; atomcnt=collections.Counter(); nterm=collections.Counter()
    kindcnt=collections.Counter()
    for idx,ln in enumerate(lines):
        ln=ln.strip()
        if not ln: continue
        e=parse_line(ln)
        m,S=core_of(e)
        if S is None: m,S=1,e
        d=spine(S)
        eqs.append((m,d))
        nterm[len(d)]+=1
        for c,a in d: atomcnt[node_str(a)]+=1
        if idx%10000==0: print(idx,time.time()-t0,file=sys.stderr)
    print('terms/core:',sorted(nterm.items())[:5], '... max',max(nterm))
    print('distinct atoms:',len(atomcnt))
    print('top atoms:',atomcnt.most_common(8))
    # atom shape census
    sh=collections.Counter()
    for a in atomcnt:
        s=re.sub(r'x\d+','X',a); s=re.sub(r'\b\d{2,}\b','C',s)
        sh[s]+=atomcnt[a]
    print('atom shapes:'); 
    for k,v in sh.most_common(30): print('   %8d  %s'%(v,k))
    pickle.dump(eqs,open(os.path.join(HERE,'eqs3.pkl'),'wb'))
    print('done',time.time()-t0)
