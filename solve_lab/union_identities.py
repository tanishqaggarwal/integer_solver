#!/usr/bin/env python3
"""CORE: the circuit obfuscates by renaming wires through identity gates (x_a - x_b
= 0) and negations (x_a + x_b = 0). Union-find ALL such 2-term equalities across the
46275 atoms to collapse forced-equal variables. Report class sizes (esp. x_15's,
x_9770's, twist vars'), and how many atoms simplify. If a huge fraction collapses,
the real system is far smaller than 38748 vars."""
import json
from collections import defaultdict
from propagate import load_atoms, atom_vars

def main():
    A = load_atoms()
    control = set(json.load(open('control_bits.json')))
    parent = {}; sign = {}   # union-find with sign (x = sign * rep)
    def find(x):
        if x not in parent: parent[x]=x; sign[x]=1
        root=x; s=1
        while parent[root]!=root:
            s*=sign[root]; root=parent[root]
        # path compress
        cur=x; sc=1
        while parent[cur]!=root:
            nxt=parent[cur]; ns=sign[cur]
            parent[cur]=root; sign[cur]=s* (1)  # recompute below
            cur=nxt
        return root, s
    # simpler union-find without full sign compression correctness edge-cases:
    parent={}; rank={}
    def f(x):
        parent.setdefault(x,x)
        while parent[x]!=x:
            parent[x]=parent[parent[x]]; x=parent[x]
        return x
    def u(a,b):
        ra,rb=f(a),f(b)
        if ra!=rb: parent[ra]=rb

    nid=0
    for poly in A:
        if len(poly)==2:
            (m1,c1),(m2,c2)=list(poly.items())
            if len(m1)==1 and len(m2)==1:   # c1*x_a + c2*x_b = 0
                a,b=m1[0],m2[0]
                if abs(c1)==abs(c2):        # x_a = +/- x_b  (identity or negation)
                    u(a,b); nid+=1
    # classes
    classes=defaultdict(list)
    for x in list(parent): classes[f(x)].append(x)
    sizes=sorted((len(v) for v in classes.values()), reverse=True)
    print(f"identity/negation 2-term atoms: {nid}")
    print(f"variables in equalities: {len(parent)}; classes: {len(classes)}")
    print(f"largest class sizes: {sizes[:15]}")
    big=[v for v in classes.values() if len(v)>=5]
    print(f"classes with >=5 members: {len(big)}; vars in them: {sum(len(v) for v in big)}")
    # x_15's class, twist vars' classes
    for t in [15,9770,3183,18274,17728,414,37917,24026,38215]:
        cls=classes.get(f(t),[t])
        print(f"  x_{t}: class size {len(cls)}; sample {sorted(cls)[:8]}")
    # how many atoms become trivial/simpler after merging reps?
    rep=lambda x: f(x)
    trivial=0; reduced_vars=set()
    for poly in A:
        vs=set(rep(v) for v in atom_vars(poly))
        reduced_vars|=vs
    print(f"distinct REP variables across all atoms: {len(reduced_vars)} (was {len(set().union(*[atom_vars(p) for p in A if p]))})")

if __name__=='__main__':
    main()
