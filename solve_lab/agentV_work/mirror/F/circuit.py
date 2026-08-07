#!/usr/bin/env python3
"""Atom-level circuit model."""
import sys, os, pickle, collections, re
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE)
from parse import node_str, const_val

def vars_of(n, s):
    o=n[0]
    if o=='v': s.add(n[1])
    elif o=='c': pass
    elif o=='neg': vars_of(n[1],s)
    else: vars_of(n[1],s); vars_of(n[2],s)
    return s

def classify(a):
    """Return dict describing atom a."""
    s=node_str(a)
    V=vars_of(a,set())
    # definition: top node is '-' with left a bare var not in right
    if a[0]=='-' and a[1][0]=='v':
        lv=a[1][1]; rv=vars_of(a[2],set())
        if lv not in rv:
            return ('def', lv, a[2], V)
    # X*(X-1) / X*(1-X) style boolean
    return ('other', None, a, V)

def main():
    eqs=pickle.load(open(os.path.join(HERE,'eqs3.pkl'),'rb'))
    atoms={}   # str -> (ast, class)
    eqatoms=[] # per eq: list of (coef, atomstr)
    for m,d in eqs:
        row=[]
        for c,a in d:
            s=node_str(a)
            if s not in atoms: atoms[s]=classify(a)
            row.append((c,s))
        eqatoms.append((m,row))
    print('atoms',len(atoms))
    defs=collections.defaultdict(list)
    others=[]
    for s,(k,lv,rhs,V) in atoms.items():
        if k=='def': defs[lv].append(s)
        else: others.append(s)
    print('def-atoms defining %d distinct vars; non-def atoms %d'%(len(defs),len(others)))
    multi={v:ss for v,ss in defs.items() if len(ss)>1}
    print('vars with >1 definition atom:',len(multi))
    # census of non-def atoms
    sh=collections.Counter()
    for s in others:
        t=re.sub(r'x\d+','X',s); t=re.sub(r'\b\d{2,}\b','C',t); sh[t]+=1
    print('non-def atom shapes:', sh.most_common(20))
    pickle.dump((atoms,eqatoms),open(os.path.join(HERE,'circ.pkl'),'wb'))
    # how many equations consist ONLY of def atoms with distinct defined vars
    return atoms,eqatoms,defs,others

if __name__=='__main__': main()
