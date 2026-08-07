#!/usr/bin/env python3
"""Strict spine + atom classification + definition graph."""
import sys,os,pickle,collections,re,time
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE)
from parse import parse_line,node_str,const_val
from parse2 import core_of,factors

def spine(S):
    terms=[]; n=S
    while n[0]=='+' and n[2][0]=='*':
        terms.append(n[2]); n=n[1]
    terms.append(n); terms.reverse(); return terms

def split_coef(nd):
    fs=[]; factors(nd,fs); k=1; nonc=[]
    for f in fs:
        cv=const_val(f)
        if cv is None: nonc.append(f)
        else: k*=cv
    if not nonc: return k,('c',1)
    if len(nonc)==1: return k,nonc[0]
    n2=nonc[0]
    for f in nonc[1:]: n2=('*',n2,f)
    return k,n2

def vars_of(n,s=None):
    if s is None: s=set()
    o=n[0]
    if o=='v': s.add(n[1])
    elif o=='c': pass
    elif o=='neg': vars_of(n[1],s)
    else: vars_of(n[1],s); vars_of(n[2],s)
    return s

def evalnode(n,v):
    o=n[0]
    if o=='v': return v[n[1]]
    if o=='c': return n[1]
    if o=='neg': return -evalnode(n[1],v)
    a=evalnode(n[1],v); b=evalnode(n[2],v)
    return a+b if o=='+' else (a-b if o=='-' else a*b)

def classify(a):
    """('def',var,rhs) | ('bool',var) | ('cons',ast)"""
    if a[0]=='-' and a[1][0]=='v':
        lv=a[1][1]
        if lv not in vars_of(a[2]): return ('def',lv,a[2])
    if a[0]=='-' and a[2][0]=='v':
        rv=a[2][1]
        if rv not in vars_of(a[1]): return ('def',rv,('neg',a[1]))  # X*X - c  => c = X*X  (sign: atom = L - rv, zero iff rv=L)
    return ('cons',None,a)

def main():
    t0=time.time()
    lines=open(os.path.join(HERE,'..','..','EQUATIONS.txt')).read().splitlines()
    atoms={}; eqrows=[]
    for i,ln in enumerate(lines):
        ln=ln.strip()
        if not ln: continue
        e=parse_line(ln); m,S=core_of(e)
        if S is None: S=e
        row=[]
        for nd in spine(S):
            k,a=split_coef(nd)
            s=node_str(a)
            if s not in atoms: atoms[s]=a
            row.append((k,s))
        eqrows.append(row)
        if i%10000==0: print(i,time.time()-t0,file=sys.stderr)
    print('n atoms',len(atoms))
    cls={s:classify(a) for s,a in atoms.items()}
    defs=collections.defaultdict(list); cons=[]
    for s,c in cls.items():
        if c[0]=='def': defs[c[1]].append(s)
        else: cons.append(s)
    print('defined vars',len(defs),'constraint atoms',len(cons))
    multi=sorted((v for v,ss in defs.items() if len(ss)>1))
    print('multi-defined vars',len(multi))
    ct=collections.Counter(len(ss) for ss in defs.values())
    print('defs-per-var histogram',sorted(ct.items()))
    sh=collections.Counter()
    for s in cons:
        t=re.sub(r'x\d+','X',s); t=re.sub(r'\d{2,}','C',t); sh[t]+=1
    print('constraint atom shapes:',sh.most_common(25))
    pickle.dump({'atoms':atoms,'cls':cls,'eqrows':eqrows},open(os.path.join(HERE,'circ2.pkl'),'wb'))
    print('done',time.time()-t0)

if __name__=='__main__': main()
