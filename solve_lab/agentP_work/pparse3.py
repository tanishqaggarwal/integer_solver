#!/usr/bin/env python3
"""Agent P independent parse v3 — keeps parenthesization so syntactic ATOMS survive."""
import re, sys, pickle
from collections import defaultdict
sys.setrecursionlimit(200000)
EQ='/home/user/integer_solver/EQUATIONS.txt'
TOK=re.compile(r'x_\d+|\d+|[()+\-*]')

def tokenize(s): return TOK.findall(s)

class P:
    def __init__(s,t): s.t=t; s.i=0
    def peek(s): return s.t[s.i] if s.i<len(s.t) else None
    def eat(s,x=None):
        v=s.t[s.i]
        if x and v!=x: raise SyntaxError(f"{x}!={v}")
        s.i+=1; return v
    def expr(s):
        n=s.term()
        while s.peek() in ('+','-'):
            o=s.eat(); r=s.term()
            n=('+',n,r) if o=='+' else ('-',n,r)
        return n
    def term(s):
        n=s.un()
        while s.peek()=='*': s.eat(); n=('*',n,s.un())
        return n
    def un(s):
        if s.peek()=='-': s.eat(); return ('neg',s.un())
        if s.peek()=='+': s.eat(); return s.un()
        return s.at()
    def at(s):
        t=s.peek()
        if t=='(':
            s.eat('('); n=s.expr(); s.eat(')')
            return n if n[0]=='p' else ('p',n)     # collapse redundant parens
        s.eat()
        return ('v',int(t[2:])) if t.startswith('x_') else ('c',int(t))

def parse(x):
    p=P(tokenize(x)); n=p.expr()
    assert p.i==len(p.t)
    return n

def unp(n):
    while n[0]=='p': n=n[1]
    return n

# --- flatten add chain WITHOUT descending into parens ---
def fadd(n,sg=1,out=None):
    if out is None: out=[]
    if n[0]=='+': fadd(n[1],sg,out); fadd(n[2],sg,out)
    elif n[0]=='-': fadd(n[1],sg,out); fadd(n[2],-sg,out)
    elif n[0]=='neg': fadd(n[1],-sg,out)
    else: out.append((sg,n))
    return out

def fmul(n,out=None):
    if out is None: out=[]
    if n[0]=='*': fmul(n[1],out); fmul(n[2],out)
    elif n[0]=='neg': out.append(('NEG',)); fmul(n[1],out)
    else: out.append(n)
    return out

def pmul(a,b):
    r=defaultdict(int)
    for m1,c1 in a.items():
        for m2,c2 in b.items(): r[tuple(sorted(m1+m2))]+=c1*c2
    return {k:x for k,x in r.items() if x}
def padd(a,b):
    r=dict(a)
    for m,c in b.items():
        r[m]=r.get(m,0)+c
        if not r[m]: del r[m]
    return r
def topoly(n):
    t=n[0]
    if t=='p': return topoly(n[1])
    if t=='c': return {():n[1]} if n[1] else {}
    if t=='v': return {(n[1],):1}
    if t=='+': return padd(topoly(n[1]),topoly(n[2]))
    if t=='-': return padd(topoly(n[1]),{m:-c for m,c in topoly(n[2]).items()})
    if t=='neg': return {m:-c for m,c in topoly(n[1]).items()}
    if t=='*': return pmul(topoly(n[1]),topoly(n[2]))
    raise ValueError(t)
def isc(n):
    n=unp(n)
    return n[0]=='c' or (n[0]=='neg' and isc(n[1]))
def cval(n):
    n=unp(n)
    return n[1] if n[0]=='c' else -cval(n[1])
def key(p): return tuple(sorted(p.items()))

def split(node):
    """summand -> (scalar, [nonconst factor nodes])"""
    sc=1; nc=[]
    for f in fmul(node):
        if f==('NEG',): sc=-sc
        elif isc(f): sc*=cval(f)
        else: nc.append(f)
    return sc,nc

def peel(ast):
    scal=1; pw=1; node=ast
    while True:
        node2=unp(node)
        parts=[]
        for sg,nd in fadd(node2):
            sc,nc=split(nd); parts.append((sg*sc,nc))
        sigs=set(tuple(sorted(key(topoly(f)) for f in nc)) for sc,nc in parts)
        if len(sigs)!=1: break
        nc0=parts[0][1]
        if not nc0:
            return scal*sum(s for s,_ in parts),0,None
        if len(set(key(topoly(f)) for f in nc0))!=1: break
        if len(parts)==1 and len(nc0)==1 and parts[0][0]==1 and unp(nc0[0]) is node2: break
        scal*=sum(s for s,_ in parts); pw*=len(nc0); node=nc0[0]
    return scal,pw,node

def main():
    lines=[l.strip() for l in open(EQ) if l.strip()]
    a2i={}; AP=[]; ATXT=[]
    rows=[]
    for ei,line in enumerate(lines):
        scal,pw,L=peel(parse(line.rsplit('=',1)[0]))
        row=[]
        if L is not None:
            for sg,nd in fadd(unp(L)):
                sc,nc=split(nd)
                ap={():1}
                for f in nc: ap=pmul(ap,topoly(f))
                k=key(ap)
                i=a2i.get(k)
                if i is None:
                    i=len(AP); a2i[k]=i; AP.append(ap); ATXT.append(nc)
                row.append((sg*sc,i))
        rows.append({'scal':scal,'pw':pw,'row':row})
    print("eqs",len(rows),"distinct atoms",len(AP))
    from collections import Counter
    print("pw:",Counter(r['pw'] for r in rows))
    print("atoms/eq:",sorted(Counter(len(r['row']) for r in rows).items()))
    u=Counter()
    for r in rows:
        for c,a in r['row']: u[a]+=1
    print("atom usage hist:",sorted(Counter(u.values()).items())[:12])
    dg=Counter(); nv=Counter()
    for ap in AP:
        dg[max((len(m) for m in ap),default=0)]+=1
        s=set()
        for m in ap: s.update(m)
        nv[len(s)]+=1
    print("atom deg:",dict(sorted(dg.items())))
    print("atom nvars:",dict(sorted(nv.items())))
    pickle.dump({'rows':rows,'AP':AP},open('/home/user/integer_solver/solve_lab/agentP_work/model3.pkl','wb'))
    print("saved model3.pkl")

if __name__=='__main__': main()
