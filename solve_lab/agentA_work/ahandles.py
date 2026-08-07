"""Exact private-handle census: for each atom, which variables occur ONLY in it,
and the lattice (gcd of linear coefficients at the current point) of shifts they give."""
import sys, json, collections; sys.path.insert(0,'.')
import env, lib as L
from math import gcd
P=env.P
v=env.load_best()

def build(v):
    solo = collections.defaultdict(list)   # atom -> [vars occurring only in this atom]
    for u,ats in L.var_atoms.items():
        if len(ats)==1: solo[ats[0]].append(u)
    gran={}
    for a,us in solo.items():
        Pp=L.polys[a]; g=0
        for u in us:
            # coefficient of u treating others at current value; must be linear in u
            c=0; ok=True
            for m,cc in Pp.items():
                if m.count(u)>1: ok=False; break
                if u in m:
                    t=cc
                    for w in m:
                        if w!=u: t*=v[w]
                    c+=t
            if ok and c: g=gcd(g,abs(c))
        gran[a]=g
    return solo,gran

if __name__=='__main__':
    solo,gran=build(v)
    print('atoms with >=1 private var: %d'%len(solo))
    hist=collections.Counter()
    for a,g in gran.items():
        hist['p' if g==P else ('0' if g==0 else ('1' if g==1 else str(g) if g<10**9 else 'other'))]+=1
    print('granularity histogram:', dict(hist))
    json.dump({str(a):[solo[a],str(gran[a])] for a in solo}, open('handles.json','w'))
    # which of the residual-region atoms have private handles?
    for a in env.SEVEN+[22231,37887,29090,39166,40066,40932,40005,40121]:
        print('  a%-6d private=%s gran=%s'%(a, solo.get(a,[]), 'p' if gran.get(a)==P else gran.get(a)))
