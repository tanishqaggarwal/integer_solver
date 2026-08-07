#!/usr/bin/env python3
import re, sys, json, pickle, time, os, collections
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0,HERE)
from parse import parse_line, const_val, node_str, flatten_sum
EQ = os.path.join(HERE, '..', '..', 'EQUATIONS.txt')

def factors(n, out):
    if n[0]=='*': factors(n[1],out); factors(n[2],out)
    elif n[0]=='neg': out.append(('c',-1)); factors(n[1],out)
    else: out.append(n)

def core_of(e):
    """Return (mult, core_node)."""
    terms=[]; flatten_sum_top(e, 1, terms)
    cores=[]
    for sg,nd in terms:
        fs=[]; factors(nd,fs)
        k=sg; nonc=[]
        for f in fs:
            cv=const_val(f)
            if cv is None: nonc.append(f)
            else: k*=cv
        cores.append((k,nonc))
    # all cores must reference the same node
    if len(cores)==1 and len(cores[0][1])==1:
        return cores[0][0], cores[0][1][0]
    if len(cores)==1 and len(cores[0][1])>=2:
        strs=set(node_str(x) for x in cores[0][1])
        if len(strs)==1:
            return None, cores[0][1][0]   # power
        return None, None
    if len(cores)>1:
        strs=set()
        for k,nc in cores:
            if len(nc)!=1: return None,None
            strs.add(node_str(nc[0]))
        if len(strs)==1:
            return sum(k for k,_ in cores), cores[0][1][0]
    return None, None

def flatten_sum_top(n, sign, out):
    if n[0]=='+': flatten_sum_top(n[1],sign,out); flatten_sum_top(n[2],sign,out)
    elif n[0]=='-': flatten_sum_top(n[1],sign,out); flatten_sum_top(n[2],-sign,out)
    elif n[0]=='neg': flatten_sum_top(n[1],-sign,out)
    else: out.append((sign,n))

def decomp_core(S):
    """S = sum_i c_i * A_i  (top-level)."""
    terms=[]; flatten_sum_top(S,1,terms)
    res=[]
    for sg,nd in terms:
        fs=[]; factors(nd,fs)
        k=sg; nonc=[]
        for f in fs:
            cv=const_val(f)
            if cv is None: nonc.append(f)
            else: k*=cv
        if len(nonc)==0: res.append((k,('c',1)))
        elif len(nonc)==1: res.append((k,nonc[0]))
        else:
            # product atom, keep as-is
            nd2=nonc[0]
            for f in nonc[1:]: nd2=('*',nd2,f)
            res.append((k,nd2))
    return res

if __name__=='__main__':
    t0=time.time()
    lines=open(EQ).read().splitlines()
    eqs=[]
    bad=0
    shapecnt=collections.Counter()
    atomcnt=collections.Counter()
    for idx,ln in enumerate(lines):
        ln=ln.strip()
        if not ln: continue
        e=parse_line(ln)
        m,S=core_of(e)
        if S is None:
            # LHS is itself the core
            m,S=1,e
            bad+=1
        d=decomp_core(S)
        eqs.append((m,d))
        shapecnt[len(d)]+=1
        for c,a in d: atomcnt[node_str(a)]+=1
        if idx%10000==0: print(idx,time.time()-t0,file=sys.stderr)
    print('no-outer-mult (core=LHS):',bad)
    print('num terms per core, top:',shapecnt.most_common(15))
    print('distinct atom strings:',len(atomcnt))
    print('top atoms:',atomcnt.most_common(10))
    pickle.dump(eqs,open(os.path.join(HERE,'eqs.pkl'),'wb'))
    print('done',time.time()-t0)
