#!/usr/bin/env python3
"""Agent P independent parse v4 — correct chain/base split, true syntactic atoms."""
import sys, pickle
from collections import Counter, defaultdict
sys.setrecursionlimit(200000)
from pparse3 import parse, unp, split, topoly, key, peel, fmul, isc, cval

EQ='/home/user/integer_solver/EQUATIONS.txt'

def is_chain_step(right):
    """right summand is `(const)*(compound gate)` -> a chain step, not part of the base."""
    n=right
    if n[0]!='*': return False
    sc,nc=split(n)
    if len(nc)!=1: return False
    c=unp(nc[0])
    if c[0]=='v' or c[0]=='c': return False      # bare var/const -> part of base
    return True

def decompose(L):
    out=[]
    n=unp(L)
    while True:
        if n[0] in ('+','-') and is_chain_step(n[2]):
            out.append((1 if n[0]=='+' else -1, n[2]))
            left=n[1]
            n=unp(left)
        else:
            out.append((1, n)); break
    out.reverse()
    return out

def main():
    lines=[l.strip() for l in open(EQ) if l.strip()]
    a2i={}; AP=[]; AN=[]
    rows=[]
    for ei,line in enumerate(lines):
        scal,pw,L=peel(parse(line.rsplit('=',1)[0]))
        row=[]
        if L is not None:
            for sg,nd in decompose(L):
                sc,nc=split(nd)
                ap={():1}
                for f in nc: ap=pmul_(ap,topoly(f))
                k=key(ap); i=a2i.get(k)
                if i is None:
                    i=len(AP); a2i[k]=i; AP.append(ap); AN.append(nc)
                row.append((sg*sc,i))
        rows.append({'scal':scal,'pw':pw,'row':row})
    print("eqs",len(rows),"distinct atoms",len(AP))
    print("atoms/eq:",sorted(Counter(len(r['row']) for r in rows).items()))
    u=Counter()
    for r in rows:
        for c,a in r['row']: u[a]+=1
    print("atom usage hist:",sorted(Counter(u.values()).items()))
    nv=Counter(); dg=Counter()
    for ap in AP:
        s=set()
        for m in ap: s.update(m)
        nv[len(s)]+=1; dg[max((len(m) for m in ap),default=0)]+=1
    print("atom nvars:",dict(sorted(nv.items())))
    print("atom deg:",dict(sorted(dg.items())))
    pickle.dump({'rows':rows,'AP':AP},open('/home/user/integer_solver/solve_lab/agentP_work/model4.pkl','wb'))
    print("saved model4.pkl")

from pparse3 import pmul as pmul_
if __name__=='__main__': main()
