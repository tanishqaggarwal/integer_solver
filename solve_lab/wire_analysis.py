#!/usr/bin/env python3
"""CORE: substitute the 220-var wire (all = +/- V = +/- x_15) into every atom and
determine whether V != 0 is consistent. Signed union-find gives each member's sign.
Then, holding all NON-wire vars at best_partial, each atom becomes a polynomial in V:
  sum over monomials -> collect by power of V.  Report atoms that force V (linear
c*V + c0 = 0 with c0 from best-partial != 0 -> V determined; conflicts => V=0), and
whether any atom is c*V=0 or V^2*... forcing V=0. If V is consistently determined to
one nonzero value, the escape is activatable."""
import json
from collections import defaultdict, Counter
from propagate import load_atoms, atom_vars

def main():
    A = load_atoms()
    best = json.load(open('best/best_partial_39019.json')); bv={int(k[2:]):v for k,v in best.items()}
    def val(v): return bv.get(v,0)
    # signed union-find: rep + sign so that x = sign * rep
    parent={}; sgn={}
    def find(x):
        parent.setdefault(x,x); sgn.setdefault(x,1)
        if parent[x]==x: return x,1
        r,s=find(parent[x])
        parent[x]=r; sgn[x]=sgn[x]*s
        return r, sgn[x]
    def union(a,b,rel):  # x_a = rel * x_b, rel in {+1,-1}
        ra,sa=find(a); rb,sb=find(b)
        if ra==rb: return
        # sa*ra_val = a_val; want a = rel*b => sa*ra = rel*sb*rb => ra = (rel*sb/sa)*rb
        parent[ra]=rb; sgn[ra]= rel*sb*sa
    for poly in A:
        if len(poly)==2:
            (m1,c1),(m2,c2)=list(poly.items())
            if len(m1)==1 and len(m2)==1 and abs(c1)==abs(c2):
                rel = -1 if (c1>0)==(c2>0) else 1   # c1*a+c2*b=0 => a = -(c2/c1) b
                union(m1[0],m2[0],rel)
    # find x_15's class + signs
    r15,_=find(15)
    wire={}
    for x in list(parent):
        r,s=find(x)
        if r==r15: wire[x]=s   # x = s * V
    print(f"wire size {len(wire)} (x=+/-V); +V:{sum(1 for s in wire.values() if s>0)}, -V:{sum(1 for s in wire.values() if s<0)}")

    # substitute wire->V, other vars at best; per atom collect poly in V
    forceV0=[]; linear=[]; determined=Counter()
    touched=0
    for a,poly in enumerate(A):
        if not (atom_vars(poly) & set(wire)): continue
        touched+=1
        # coefficient of V^k
        coefs=defaultdict(int)
        for m,c in poly.items():
            k=0; t=c
            for x in m:
                if x in wire: k+=1; t*=wire[x]
                else: t*=val(x)
            coefs[k]+=t
        c0=coefs.get(0,0); c1=coefs.get(1,0); hi=sum(1 for k in coefs if k>=2 and coefs[k])
        if c0==0 and c1==0 and hi==0: continue  # trivially 0
        if hi==0:
            # c1*V + c0 = 0
            if c1==0:
                if c0!=0: pass  # atom broken at V-independent -> best_partial residual (non-wire), ignore
            else:
                if c0%c1==0: determined[(-c0)//c1]+=1; linear.append((a,(-c0)//c1))
                else: linear.append((a,'noninteger'))
        else:
            # has V^2+ terms; check if forces V=0 (e.g. c*V^2=0)
            if c0==0 and c1==0:
                forceV0.append(a)
    print(f"atoms touching the wire: {touched}")
    print(f"atoms forcing V=0 (pure V^2.. =0, no const/linear): {len(forceV0)}: {forceV0[:10]}")
    print(f"atoms LINEAR in V determining a value: {len(linear)}")
    # do the linear atoms agree on V?
    vals=Counter(v for _,v in linear if v!='noninteger')
    print(f"  distinct V-values demanded by linear atoms: {len(vals)}")
    for vv,ct in vals.most_common(8):
        print(f"    V={vv}: demanded by {ct} atoms")
    noninteg=sum(1 for _,v in linear if v=='noninteger')
    print(f"  linear atoms with non-integer V demand: {noninteg}")

if __name__=='__main__':
    main()
